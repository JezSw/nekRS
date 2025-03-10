#include <cpptrace/cpptrace.hpp>

#include <atomic>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <new>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "symbols/symbols.hpp"
#include "unwind/unwind.hpp"
#include "demangle/demangle.hpp"
#include "utils/exception_type.hpp"
#include "utils/common.hpp"
#include "utils/utils.hpp"
#include "binary/object.hpp"
#include "binary/safe_dl.hpp"

#define ESC     "\033["
#define RESET   ESC "0m"
#define RED     ESC "31m"
#define GREEN   ESC "32m"
#define YELLOW  ESC "33m"
#define BLUE    ESC "34m"
#define MAGENTA ESC "35m"
#define CYAN    ESC "36m"

namespace cpptrace {
    CPPTRACE_FORCE_NO_INLINE
    raw_trace raw_trace::current(std::size_t skip) {
        return generate_raw_trace(skip + 1);
    }

    CPPTRACE_FORCE_NO_INLINE
    raw_trace raw_trace::current(std::size_t skip, std::size_t max_depth) {
        return generate_raw_trace(skip + 1, max_depth);
    }

    object_trace raw_trace::resolve_object_trace() const {
        try {
            return object_trace{detail::get_frames_object_info(frames)};
        } catch(...) { // NOSONAR
            if(!detail::should_absorb_trace_exceptions()) {
                throw;
            }
            return object_trace{};
        }
    }

    stacktrace raw_trace::resolve() const {
        try {
            std::vector<stacktrace_frame> trace = detail::resolve_frames(frames);
            for(auto& frame : trace) {
                frame.symbol = detail::demangle(frame.symbol);
            }
            return {std::move(trace)};
        } catch(...) { // NOSONAR
            if(!detail::should_absorb_trace_exceptions()) {
                throw;
            }
            return stacktrace{};
        }
    }

    void raw_trace::clear() {
        frames.clear();
    }

    bool raw_trace::empty() const noexcept {
        return frames.empty();
    }

    CPPTRACE_FORCE_NO_INLINE
    object_trace object_trace::current(std::size_t skip) {
        return generate_object_trace(skip + 1);
    }

    CPPTRACE_FORCE_NO_INLINE
    object_trace object_trace::current(std::size_t skip, std::size_t max_depth) {
        return generate_object_trace(skip + 1, max_depth);
    }

    stacktrace object_trace::resolve() const {
        try {
            std::vector<stacktrace_frame> trace = detail::resolve_frames(frames);
            for(auto& frame : trace) {
                frame.symbol = detail::demangle(frame.symbol);
            }
            return {std::move(trace)};
        } catch(...) { // NOSONAR
            if(!detail::should_absorb_trace_exceptions()) {
                throw;
            }
            return stacktrace();
        }
    }

    void object_trace::clear() {
        frames.clear();
    }

    bool object_trace::empty() const noexcept {
        return frames.empty();
    }

    std::string stacktrace_frame::to_string() const {
        std::ostringstream oss;
        oss << *this;
        return std::move(oss).str();
    }

    std::ostream& operator<<(std::ostream& stream, const stacktrace_frame& frame) {
        stream
            << std::hex
            << "0x"
            << std::setw(2 * sizeof(frame_ptr))
            << std::setfill('0')
            << frame.raw_address
            << std::dec
            << std::setfill(' ');
        if(!frame.symbol.empty()) {
            stream
                << " in "
                << frame.symbol;
        }
        if(!frame.filename.empty()) {
            stream
                << " at "
                << frame.filename;
            if(frame.line.has_value()) {
                stream
                    << ":"
                    << frame.line.value();
                if(frame.column.has_value()) {
                    stream << frame.column.value();
                }
            }
        }
        return stream;
    }

    CPPTRACE_FORCE_NO_INLINE
    stacktrace stacktrace::current(std::size_t skip) {
        return generate_trace(skip + 1);
    }

    CPPTRACE_FORCE_NO_INLINE
    stacktrace stacktrace::current(std::size_t skip, std::size_t max_depth) {
        return generate_trace(skip + 1, max_depth);
    }

    void stacktrace::print() const {
        print(std::cerr, true);
    }

    void stacktrace::print(std::ostream& stream) const {
        print(stream, true);
    }

    void stacktrace::print(std::ostream& stream, bool color) const {
        print(stream, color, true, nullptr);
    }

    void stacktrace::print(std::ostream& stream, bool color, bool newline_at_end, const char* header) const {
        if(
            color && (
                (&stream == &std::cout && isatty(stdout_fileno)) || (&stream == &std::cerr && isatty(stderr_fileno))
            )
        ) {
            detail::enable_virtual_terminal_processing_if_needed();
        }
        stream<<(header ? header : "Stack trace (most recent call first):") << '\n';
        std::size_t counter = 0;
        if(frames.empty()) {
            stream<<"<empty trace>" << '\n';
            return;
        }
        const auto reset   = color ? ESC "0m" : "";
        const auto green   = color ? ESC "32m" : "";
        const auto yellow  = color ? ESC "33m" : "";
        const auto blue    = color ? ESC "34m" : "";
        const auto frame_number_width = detail::n_digits(static_cast<int>(frames.size()) - 1);
        for(const auto& frame : frames) {
            stream
                << '#'
                << std::setw(static_cast<int>(frame_number_width))
                << std::left
                << counter
                << std::right
                << " ";
            if(frame.is_inline) {
                stream
                    << std::setw(2 * sizeof(frame_ptr) + 2)
                    << "(inlined)";
            } else {
                stream
                    << std::hex
                    << blue
                    << "0x"
                    << std::setw(2 * sizeof(frame_ptr))
                    << std::setfill('0')
                    << frame.raw_address
                    << std::dec
                    << std::setfill(' ')
                    << reset;
            }
            if(!frame.symbol.empty()) {
                stream
                    << " in "
                    << yellow
                    << frame.symbol
                    << reset;
            }
            if(!frame.filename.empty()) {
                stream
                    << " at "
                    << green
                    << frame.filename
                    << reset;
                if(frame.line.has_value()) {
                    stream
                        << ":"
                        << blue
                        << frame.line.value()
                        << reset;
                    if(frame.column.has_value()) {
                        stream << ':'
                            << blue
                            << std::to_string(frame.column.value())
                            << reset;
                    }
                }
            }
            if(newline_at_end || &frame != &frames.back()) {
                stream << '\n';
            }
            counter++;
        }
    }

    void stacktrace::clear() {
        frames.clear();
    }

    bool stacktrace::empty() const noexcept {
        return frames.empty();
    }

    std::string stacktrace::to_string(bool color) const {
        std::ostringstream oss;
        print(oss, color, false, nullptr);
        return std::move(oss).str();
    }

    std::ostream& operator<<(std::ostream& stream, const stacktrace& trace) {
        return stream << trace.to_string();
    }

    CPPTRACE_FORCE_NO_INLINE
    raw_trace generate_raw_trace(std::size_t skip) {
        try {
            return raw_trace{detail::capture_frames(skip + 1, SIZE_MAX)};
        } catch(...) { // NOSONAR
            if(!detail::should_absorb_trace_exceptions()) {
                throw;
            }
            return raw_trace{};
        }
    }

    CPPTRACE_FORCE_NO_INLINE
    raw_trace generate_raw_trace(std::size_t skip, std::size_t max_depth) {
        try {
            return raw_trace{detail::capture_frames(skip + 1, max_depth)};
        } catch(...) { // NOSONAR
            if(!detail::should_absorb_trace_exceptions()) {
                throw;
            }
            return raw_trace{};
        }
    }

    CPPTRACE_FORCE_NO_INLINE
    std::size_t safe_generate_raw_trace(frame_ptr* buffer, std::size_t size, std::size_t skip) {
        return detail::safe_capture_frames(buffer, size, skip + 1, SIZE_MAX);
    }

    CPPTRACE_FORCE_NO_INLINE
    std::size_t safe_generate_raw_trace(
         frame_ptr* buffer,
         std::size_t size,
         std::size_t skip,
         std::size_t max_depth
    ) {
        return detail::safe_capture_frames(buffer, size, skip + 1, max_depth);
    }

    CPPTRACE_FORCE_NO_INLINE
    object_trace generate_object_trace(std::size_t skip) {
        try {
            return object_trace{detail::get_frames_object_info(detail::capture_frames(skip + 1, SIZE_MAX))};
        } catch(...) { // NOSONAR
            if(!detail::should_absorb_trace_exceptions()) {
                throw;
            }
            return object_trace{};
        }
    }

    CPPTRACE_FORCE_NO_INLINE
    object_trace generate_object_trace(std::size_t skip, std::size_t max_depth) {
        try {
            return object_trace{detail::get_frames_object_info(detail::capture_frames(skip + 1, max_depth))};
        } catch(...) { // NOSONAR
            if(!detail::should_absorb_trace_exceptions()) {
                throw;
            }
            return object_trace{};
        }
    }

    CPPTRACE_FORCE_NO_INLINE
    stacktrace generate_trace(std::size_t skip) {
        return generate_trace(skip + 1, SIZE_MAX);
    }

    CPPTRACE_FORCE_NO_INLINE
    stacktrace generate_trace(std::size_t skip, std::size_t max_depth) {
        try {
            std::vector<frame_ptr> frames = detail::capture_frames(skip + 1, max_depth);
            std::vector<stacktrace_frame> trace = detail::resolve_frames(frames);
            for(auto& frame : trace) {
                frame.symbol = detail::demangle(frame.symbol);
            }
            return {std::move(trace)};
        } catch(...) { // NOSONAR
            if(!detail::should_absorb_trace_exceptions()) {
                throw;
            }
            return stacktrace();
        }
    }

    object_frame safe_object_frame::resolve() const {
        return detail::resolve_safe_object_frame(*this);
    }

    void get_safe_object_frame(frame_ptr address, safe_object_frame* out) {
        detail::get_safe_object_frame(address, out);
    }

    std::string demangle(const std::string& name) {
        return detail::demangle(name);
    }

    bool isatty(int fd) {
        return detail::isatty(fd);
    }

    extern const int stdin_fileno = detail::fileno(stdin);
    extern const int stdout_fileno = detail::fileno(stdout);
    extern const int stderr_fileno = detail::fileno(stderr);

    CPPTRACE_FORCE_NO_INLINE void print_terminate_trace() {
        generate_trace(1).print(
            std::cerr,
            isatty(stderr_fileno),
            true,
            "Stack trace to reach terminate handler (most recent call first):"
        );
    }

    [[noreturn]] void terminate_handler() {
        // TODO: Support std::nested_exception?
        try {
            auto ptr = std::current_exception();
            if(ptr == nullptr) {
                std::cerr << "terminate called without an active exception\n";
                print_terminate_trace();
            } else {
                std::rethrow_exception(ptr);
            }
        } catch(cpptrace::exception& e) {
            std::cerr << "Terminate called after throwing an instance of "
                      << demangle(typeid(e).name())
                      << ": "
                      << e.message()
                      << '\n';
            e.trace().print(std::cerr, isatty(stderr_fileno));
        } catch(std::exception& e) {
            std::cerr << "Terminate called after throwing an instance of "
                      << demangle(typeid(e).name())
                      << ": "
                      << e.what()
                      << '\n';
            print_terminate_trace();
        } catch(...) {
            std::cerr << "Terminate called after throwing an instance of "
                      << detail::exception_type_name()
                      << "\n";
            print_terminate_trace();
        }
        std::flush(std::cerr);
        abort();
    }

    void register_terminate_handler() {
        std::set_terminate(terminate_handler);
    }

    namespace detail {
        std::atomic_bool absorb_trace_exceptions(true); // NOSONAR
        std::atomic_bool resolve_inlined_calls(true); // NOSONAR
        std::atomic<enum cache_mode> cache_mode(cache_mode::prioritize_speed); // NOSONAR
    }

    void absorb_trace_exceptions(bool absorb) {
        detail::absorb_trace_exceptions = absorb;
    }

     void enable_inlined_call_resolution(bool enable) {
        detail::resolve_inlined_calls = enable;
    }

    namespace experimental {
        void set_cache_mode(cache_mode mode) {
            detail::cache_mode = mode;
        }
    }

    namespace detail {
        bool should_absorb_trace_exceptions() {
            return absorb_trace_exceptions;
        }

        bool should_resolve_inlined_calls() {
            return resolve_inlined_calls;
        }

        enum cache_mode get_cache_mode() {
            return cache_mode;
        }

        CPPTRACE_FORCE_NO_INLINE
        raw_trace get_raw_trace_and_absorb(std::size_t skip, std::size_t max_depth) noexcept {
            try {
                return generate_raw_trace(skip + 1, max_depth);
            } catch(const std::exception& e) {
                if(!detail::should_absorb_trace_exceptions()) {
                    // TODO: Append to message somehow
                    std::fprintf(
                        stderr,
                        "Cpptrace: Exception occurred while resolving trace in cpptrace::exception object:\n%s\n",
                        e.what()
                    );
                }
                return raw_trace{};
            }
        }

        CPPTRACE_FORCE_NO_INLINE
        raw_trace get_raw_trace_and_absorb(std::size_t skip) noexcept {
            return get_raw_trace_and_absorb(skip + 1, SIZE_MAX);
        }

        lazy_trace_holder::lazy_trace_holder(const lazy_trace_holder& other) : resolved(other.resolved) {
            if(other.resolved) {
                new (&resolved_trace) stacktrace(other.resolved_trace);
            } else {
                new (&trace) raw_trace(other.trace);
            }
        }
        lazy_trace_holder::lazy_trace_holder(lazy_trace_holder&& other) noexcept : resolved(other.resolved) {
            if(other.resolved) {
                new (&resolved_trace) stacktrace(std::move(other.resolved_trace));
            } else {
                new (&trace) raw_trace(std::move(other.trace));
            }
        }
        lazy_trace_holder& lazy_trace_holder::operator=(const lazy_trace_holder& other) {
            clear();
            resolved = other.resolved;
            if(other.resolved) {
                new (&resolved_trace) stacktrace(other.resolved_trace);
            } else {
                new (&trace) raw_trace(other.trace);
            }
            return *this;
        }
        lazy_trace_holder& lazy_trace_holder::operator=(lazy_trace_holder&& other) noexcept {
            clear();
            resolved = other.resolved;
            if(other.resolved) {
                new (&resolved_trace) stacktrace(std::move(other.resolved_trace));
            } else {
                new (&trace) raw_trace(std::move(other.trace));
            }
            return *this;
        }
        lazy_trace_holder::~lazy_trace_holder() {
            clear();
        }
        // access
        stacktrace& lazy_trace_holder::get_resolved_trace() {
            if(!resolved) {
                stacktrace new_trace;
                try {
                    if(resolved_trace.empty() && !trace.empty()) {
                        resolved_trace = trace.resolve();
                        trace.clear();
                    }
                    new_trace = trace.resolve();
                } catch(const std::exception& e) {
                    if(!detail::should_absorb_trace_exceptions()) {
                        // TODO: Append to message somehow?
                        std::fprintf(
                            stderr,
                            "Exception occurred while resolving trace in cpptrace::exception object:\n%s\n",
                            e.what()
                        );
                    }
                }
                trace.~raw_trace();
                new (&resolved_trace) stacktrace(std::move(new_trace));
                resolved = true;
            }
            return resolved_trace;
        }
        const stacktrace& lazy_trace_holder::get_resolved_trace() const {
            if(!resolved) {
                throw std::logic_error(
                    "cpptrace::detaillazy_trace_holder::get_resolved_trace called on unresolved const object"
                );
            }
            return resolved_trace;
        }
        void lazy_trace_holder::clear() {
            if(resolved) {
                resolved_trace.~stacktrace();
            } else {
                trace.~raw_trace();
            }
        }
    }

    const char* lazy_exception::what() const noexcept {
        if(what_string.empty()) {
            what_string = message() + std::string(":\n") + trace_holder.get_resolved_trace().to_string();
        }
        return what_string.c_str();
    }

    const char* lazy_exception::message() const noexcept {
        return "cpptrace::lazy_exception";
    }

    const stacktrace& lazy_exception::trace() const noexcept {
        return trace_holder.get_resolved_trace();
    }

    const char* exception_with_message::message() const noexcept {
        return user_message.c_str();
    }

    const char* nested_exception::message() const noexcept {
        if(message_value.empty()) {
            try {
                std::rethrow_exception(ptr);
            } catch(std::exception& e) {
                message_value = std::string("Nested exception: ") + e.what();
            } catch(...) {
                message_value = "Nested exception holding instance of " + detail::exception_type_name();
            }
        }
        return message_value.c_str();
    }

    std::exception_ptr nested_exception::nested_ptr() const noexcept {
        return ptr;
    }

    CPPTRACE_FORCE_NO_INLINE
    void rethrow_and_wrap_if_needed(std::size_t skip) {
        try {
            std::rethrow_exception(std::current_exception());
        } catch(cpptrace::exception&) {
            throw; // already a cpptrace::exception
        } catch(...) {
            throw nested_exception(std::current_exception(), detail::get_raw_trace_and_absorb(skip + 1));
        }
    }
}
