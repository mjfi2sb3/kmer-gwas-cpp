#ifndef mmap_io_h
#define mmap_io_h

#include <algorithm>
#include <iostream>
#include <cstring>
#include <vector>

// for mmap:
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>

//for converting [const char*] to [istream]
#include <istream>
#include <streambuf>
#include <string>

void get_reads(std::string fname, std::vector<std::pair<std::string, std::string>>& test_set, size_t max_size = 0);
template <typename T>
void read_vector(std::string fname, std::vector<T>& v);
template <typename T>
void write_vector(std::string fname, const std::vector<T>& v);
void handle_error(std::string msg);


template <typename T>
void write_vector(std::string fname, const std::vector<T>& v) {
    const size_t size_ = v.size()*sizeof(T);
    mode_t mode = S_IRUSR | S_IWUSR;
    int fd = open(fname.c_str(), O_RDWR | O_CREAT, mode);

    if (fd == -1)
        handle_error("open");

    if (posix_fallocate(fd, 0, size_) != 0) {
        close(fd);
        handle_error("posix_fallocate");
    }

    T* addr = static_cast<T*>(mmap(NULL, size_, PROT_WRITE, MAP_SHARED, fd, 0u));

    if (addr == MAP_FAILED) {
        close(fd);
        handle_error("mmap");
    }

    // Bulk copy instead of element-by-element
    std::memcpy(addr, v.data(), size_);

    msync(addr, size_, MS_SYNC);
    munmap(addr, size_);
    close(fd);

}


template <typename T>
void read_vector(std::string fname, std::vector<T>& v) {

    int fd = open(fname.c_str(), O_RDONLY);
    if (fd == -1)
        handle_error("open");

    struct stat sb;
    if (fstat(fd, &sb) == -1) {
        close(fd);
        handle_error("fstat");
    }

    auto length = sb.st_size;
    //std::cout << "file size: " << length << std::endl;

    T* addr = static_cast<T*>(mmap(NULL, length, PROT_READ, MAP_SHARED, fd, 0u));
    if (addr == MAP_FAILED) {
        close(fd);
        handle_error("mmap");
    }

    size_t N = length/sizeof(T);

    //std::cout << "N : " << N << std::endl;
    // Direct assignment instead of push_back loop
    v.assign(addr, addr + N);

    munmap(addr, length);
    close(fd);

}





#endif



