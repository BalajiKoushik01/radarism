#pragma once
#include <cstddef>
#include <cstdlib>
#include <new>

namespace radar {

// 64-byte aligned allocator: satisfies VOLK aligned kernels, FFTW SIMD paths
// and cache-line alignment.
template <class T, std::size_t Alignment = 64>
struct AlignedAllocator {
    using value_type = T;

    template <class U>
    struct rebind { using other = AlignedAllocator<U, Alignment>; };

    AlignedAllocator() noexcept = default;
    template <class U>
    AlignedAllocator(const AlignedAllocator<U, Alignment>&) noexcept {}

    [[nodiscard]] T* allocate(std::size_t n) {
        if (n > static_cast<std::size_t>(-1) / sizeof(T)) throw std::bad_array_new_length();
        std::size_t bytes = n * sizeof(T);
        if (bytes == 0) bytes = Alignment;
        bytes = ((bytes + Alignment - 1) / Alignment) * Alignment;  // aligned_alloc needs a multiple
        void* p = std::aligned_alloc(Alignment, bytes);
        if (!p) throw std::bad_alloc();
        return static_cast<T*>(p);
    }

    void deallocate(T* p, std::size_t) noexcept { std::free(p); }

    template <class U>
    bool operator==(const AlignedAllocator<U, Alignment>&) const noexcept { return true; }
};

}  // namespace radar
