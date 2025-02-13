#include <occa/internal/utils/string.hpp>
#include <occa/internal/lang/modes/dpcpp.hpp>
#include <occa/internal/lang/modes/okl.hpp>
#include <occa/internal/lang/modes/oklForStatement.hpp>
#include <occa/internal/lang/builtins/attributes.hpp>
#include <occa/internal/lang/builtins/types.hpp>
#include <occa/internal/lang/expr.hpp>
#include <occa/internal/lang/attribute.hpp>

namespace {

class dpcppLambda_t : public occa::lang::lambda_t {
public:
  int simd_length{-1};

  dpcppLambda_t(occa::lang::capture_t capture_, int simd_length_)
   : lambda_t(capture_), simd_length(simd_length_) {}

  dpcppLambda_t(const dpcppLambda_t& other) 
    : lambda_t(other), simd_length(other.simd_length) {}

  ~dpcppLambda_t() = default;

  bool equals(const type_t &other) const override {
    const dpcppLambda_t &other_ = other.to<dpcppLambda_t>();
    if (simd_length != other_.simd_length) return false;
    return lambda_t::equals(other);
  }

  void printDeclaration(occa::lang::printer &pout) const override {
    pout << "[";

    switch (this->capture) {
    case occa::lang::capture_t::byValue:
      pout << "=";
      break;
    case occa::lang::capture_t::byReference:
      pout << "&";
      break;
    default:
      pout << "???";
      break;
    }

    pout << "](";
    
    if (!args.empty()) {
      const std::string argIndent = pout.indentFromNewline();
      args[0]->printDeclaration(pout);
      for (std::size_t i = 1; i < args.size(); ++i) {
        pout << ",\n" << argIndent;
        args[i]->printDeclaration(pout);
      }
    }
    pout << ") ";
    
    if (0 < simd_length) {
      pout << "[[intel::reqd_sub_group_size("; 
      pout.print(simd_length);
      pout << ")]]";
    }

    pout << " {";

    pout.printNewline();
    pout.pushInlined(false);
    pout.addIndentation();

    body->print(pout);

    pout.removeIndentation();
    pout.popInlined();
    pout.printNewline();
    pout.printIndentation();
    pout << "}\n";
  }
};

}

namespace occa
{
  namespace lang
  {
    namespace okl
    {
      dpcppParser::dpcppParser(const occa::json &settings_)
          : withLauncher(settings_),
            kernel(externC),
            device("SYCL_EXTERNAL", qualifierType::custom),
            shared("auto", qualifierType::custom)
      {
        okl::addOklAttributes(*this);
        simd_length_default = settings_.get("simd_length",-1);
      }

      void dpcppParser::onClear()
      {
        launcherClear();
      }

      void dpcppParser::beforePreprocessing()
      {
        preprocessor.addCompilerDefine("OCCA_USING_GPU", "1");
      }

      void dpcppParser::beforeKernelSplit()
      {
        // addExtensions();

        // if(!success)
        //   return;
        // updateConstToConstant();

        if (!success)
          return;
        setFunctionQualifiers();

        if (!success)
          return;
        setSharedQualifiers();
      }

      void dpcppParser::afterKernelSplit()
      {
        addBarriers();

        if (!success)
          return;
        setupHeaders();

        if (!success)
          return;
        setupAtomics();

        if (!success)
          return;
        // Do this last!
        setupKernels();
      }

      std::string dpcppParser::getOuterIterator(const int loopIndex)
      {
        return "item_.get_group(" + occa::toString(dpcppDimensionOrder(loopIndex)) + ")";
      }

      std::string dpcppParser::getInnerIterator(const int loopIndex)
      {
        return "item_.get_local_id(" + occa::toString(dpcppDimensionOrder(loopIndex)) + ")";
      }

      std::string dpcppParser::launchBoundsAttribute(const int innerDims[3])
      {
        return "";
      }

      void dpcppParser::setupHeaders()
      {
        root.addFirst(
            *(new directiveStatement(
                &root,
                directiveToken(root.source->origin, "include <sycl/sycl.hpp>\n using namespace sycl;\n"))));
      }

      void dpcppParser::addExtensions()
      {
        if (!settings.has("extensions"))
        {
          return;
        }

        occa::json &extensions = settings["extensions"];
        if (!extensions.isObject())
        {
          return;
        }

        // @todo: Enable dpcpp extensions

        // jsonObject &extensionObj = extensions.object();
        // jsonObject::iterator it = extensionObj.begin();
        // while (it != extensionObj.end()) {
        //   const std::string &extension = it->first;
        //   const bool enabled = it->second;
        //   if (enabled) {
        //     root.addFirst(
        //       *(new pragmaStatement(
        //           &root,
        //           pragmaToken(root.source->origin,
        //                       "OPENCL EXTENSION "+ extension + " : enable\n")
        //         ))
        //     );
        //   }
        //   ++it;
        // }
      }

      void dpcppParser::addBarriers()
      {
        statementArray::from(root)
            .flatFilterByStatementType(statementType::empty, "barrier")
            .forEach([&](statement_t *smnt)
                     {
                       emptyStatement &emptySmnt = (emptyStatement &)*smnt;

                       statement_t &barrierSmnt = (*(new sourceCodeStatement(
                           emptySmnt.up,
                           emptySmnt.source,
                           "item_.barrier(sycl::access::fence_space::local_space);")));

                       emptySmnt.replaceWith(barrierSmnt);

                       delete &emptySmnt;
                     });
      }

      void dpcppParser::setupKernels()
      {
        root.children
            .filterByStatementType(
                statementType::functionDecl | statementType::function,
                "kernel")
            .forEach([&](statement_t *smnt)
                     {
                       function_t *function;

                       if (smnt->type() & statementType::functionDecl)
                       {
                         functionDeclStatement &k = ((functionDeclStatement &)*smnt);
                         function = &(k.function());

                         migrateLocalDecls(k);
                         if (!success) return;

                         int simd_length = simd_length_default;
                         if (k.hasAttribute("simd_length")) {
                           const attributeToken_t& attr = k.attributes["simd_length"];
                           simd_length = attr.args[0].expr->evaluate();
                         }

                         variable_t sycl_nditem(syclNdItem, "item_");

                         dpcppLambda_t& sycl_kernel = *(new dpcppLambda_t(capture_t::byValue, simd_length));
                         sycl_kernel.addArgument(sycl_nditem);
                         sycl_kernel.body->swap(k);

                         lambdaNode sycl_kernel_node(sycl_kernel.source, sycl_kernel);

                         variable_t sycl_ndrange(syclNdRange, "range_");
                         sycl_ndrange += pointer_t();

                         leftUnaryOpNode sycl_ndrange_node(
                             sycl_ndrange.source,
                             op::dereference,
                             variableNode(sycl_ndrange.source, sycl_ndrange.clone()));

                         exprNodeVector parallelfor_args;
                         parallelfor_args.push_back(&sycl_ndrange_node);
                         parallelfor_args.push_back(&sycl_kernel_node);

                         std::string parallelfor_name = "parallel_for<class ";
                         std::string_view function_name = function->name();
                         function_name.remove_prefix(6); // Remove _occa_ from start of kernel name
                         parallelfor_name += function_name;

                         const std::string hash_string = settings.get<std::string>("hash");
                         if (!hash_string.empty()) {
                          parallelfor_name += hash_string;
                         }
                         parallelfor_name += ">";

                         identifierNode parallelfor_node(
                             new identifierToken(originSource::builtin, "parfor"),
                             parallelfor_name);

                         callNode parallelfor_call_node(
                             parallelfor_node.token,
                             parallelfor_node,
                             parallelfor_args);

                         variable_t sycl_queue(syclQueue, "queue_");
                         sycl_queue += pointer_t();

                         binaryOpNode q_parallelfor(
                             sycl_queue.source,
                             op::arrow,
                             variableNode(sycl_queue.source, sycl_queue.clone()),
                             parallelfor_call_node);

                         k.addFirst(*(new expressionStatement(nullptr, q_parallelfor)));

                         function->addArgumentFirst(sycl_ndrange);
                         function->addArgumentFirst(sycl_queue);
                       }
                       else
                       {
                         function = &(((functionStatement *)smnt)->function());
                       }
                       setKernelQualifiers(*function);
                     });
      }

      void dpcppParser::setFunctionQualifiers()
      {
        root.children
            .filterByStatementType(statementType::functionDecl)
            .forEach([&](statement_t *smnt)
                     {
                       functionDeclStatement &funcDeclSmnt = (functionDeclStatement &)*smnt;
                       if (funcDeclSmnt.hasAttribute("kernel")) return;

                       vartype_t &vartype = funcDeclSmnt.function().returnType;
                       if (vartype.has(occa::lang::static_)) return;
                       
                       // Only add SYCL_EXTERNAL if we have external linkage
                       vartype.qualifiers.addFirst(vartype.origin(), device);
                     });
      }

      void dpcppParser::setSharedQualifiers()
      {
        statementArray::from(root).nestedForEachDeclaration(
            [&](variableDeclaration &decl, declarationStatement &declSmnt)
            {
              variable_t &var = decl.variable();

              if (var.hasAttribute("shared"))
              {
                auto *shared_value = new dpcppLocalMemoryNode(var.source->clone(),
                                                              var.vartype,
                                                              "item_");

                decl.setValue(shared_value);
                var.vartype.setType(auto_);
                var.vartype.setReferenceToken(var.source);
                var.vartype.arrays.clear();
                var.vartype.qualifiers.clear();
              }
            });
      }

      void dpcppParser::setKernelQualifiers(function_t &function)
      {
        function.returnType.add(0, kernel);

        for (auto arg : function.args)
        {
          vartype_t &type = arg->vartype;
          type = type.flatten();
          if (!(type.isPointerType() || type.referenceToken))
          {
            type.setReferenceToken(arg->source);
          }
        }
      }

      void dpcppParser::migrateLocalDecls(functionDeclStatement &kernelSmnt)
      {
        statementArray::from(kernelSmnt)
            .nestedForEachDeclaration([&](variableDeclaration &decl, declarationStatement &declSmnt)
                                      {
                                        variable_t &var = decl.variable();
                                        if (var.hasAttribute("shared"))
                                        {
                                          declSmnt.removeFromParent();
                                          kernelSmnt.addFirst(declSmnt);
                                        }
                                      });
      }

      void dpcppParser::setupAtomics()
      {
        success &= attributes::atomic::applyCodeTransformation(
            root,
            transformAtomicBlockStatement,
            transformAtomicBasicExpressionStatement);
      }

      bool dpcppParser::transformAtomicBlockStatement(blockStatement &blockSmnt)
      {
        bool transform_successful{true};
        statementArray::from(blockSmnt)
            .flatFilterByStatementType(statementType::expression)
            .forEach([&](statement_t *smnt)
                     {
                       expressionStatement &exprSmnt = static_cast<expressionStatement &>(*smnt);

                       // if(!transformAtomicBasicExpressionStatement(exprSmnt))
                       // {
                       //   transform_successful = false;
                       //   return;
                       // }
                       transformAtomicBasicExpressionStatement(exprSmnt);
                     });

        return transform_successful;
      }

      bool dpcppParser::transformAtomicBasicExpressionStatement(expressionStatement &exprSmnt)
      {
        expressionStatement &atomicSmnt = dynamic_cast<expressionStatement &>(exprSmnt.clone());

        const opType_t &opType = expr(atomicSmnt.expr).opType();

        exprNode *variable_node{nullptr};
        if (opType & operatorType::unary)
        {
          if (opType & operatorType::leftUnary)
          {
            variable_node = ((leftUnaryOpNode *)atomicSmnt.expr)->value;
          }
          else if (opType & operatorType::rightUnary)
          {
            variable_node = ((rightUnaryOpNode *)atomicSmnt.expr)->value;
          }
        }
        else if (opType & operatorType::binary)
        {
          binaryOpNode &binaryNode = *static_cast<binaryOpNode *>(atomicSmnt.expr);
          variable_node = binaryNode.leftValue;
        }
        else
        {
          atomicSmnt.printError("Unable to transform @atomic code");
          return false;
        }

        variable_t &atomic_var = *(variable_node->getVariable());
        vartype_t atomic_type = atomic_var.vartype;

        auto *atomic_ref = new dpcppAtomicNode(atomic_var.source, atomic_type, *variable_node);
        atomicSmnt.replaceExprNode(variable_node, atomic_ref);

        exprSmnt.replaceWith(atomicSmnt);
        delete &exprSmnt;

        return true;
      }

    } // namespace okl
  }   // namespace lang
} // namespace occa
