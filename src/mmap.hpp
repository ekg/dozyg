#include <fcntl.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <iostream>

namespace gyeet {

// memory mapped buffer struct
template<typename T>
struct mmap_buffer_t {
    int fd;
    off_t size;
    T *data;
    std::vector<T>::iterator begin(void) {
        return (std::vector<T>::iterator)(data + size);
    }
};

// utilities used by mmmultimap
template<typename T>
mmap_buffer_t<T> open_mmap_buffer(const char* path) {
    mmap_buffer_t<T> buffer;
    buffer.data = nullptr;
    buffer.fd = open(path, O_RDWR);
    if (buffer.fd == -1) {
        std::cerr << "[gyeet/mmap.hpp] Error: unable to open file " << path << std::endl;
        exit(1);
    }
    struct stat stats;
    if (-1 == fstat(buffer.fd, &stats)) {
        std::cerr << "[gyeet/mmap.hpp] Error: unable to fstat file " << path << std::endl;
        exit(1);
    }
    if (!(buffer.data = (T*)mmap(nullptr,
                                 stats.st_size,
                                 PROT_READ | PROT_WRITE,
                                 MAP_SHARED,
                                 buffer.fd,
                                 0
              ))) {
        std::cerr << "[gyeet/mmap.hpp] Error: unable to map buffer to file " << path << std::endl;
        exit(1);
    }
    //madvise(buffer, stats.st_size, POSIX_MADV_WILLNEED | POSIX_MADV_SEQUENTIAL);
    buffer.size = stats.st_size;
    return buffer;
}

template<typename T>
void close_mmap_buffer(const mmap_buffer_t<T>& buffer) {
    if (buffer.data != nullptr) {
        munmap((void*)buffer.data, buffer.size);
        buffer.data = nullptr;
        buffer.size = 0;
    }
    if (buffer.fd) {
        close(buffer.fd);
        buffer.fd = 0;
    }
}

}
