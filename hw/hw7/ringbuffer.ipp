#include "ringbuffer.hpp"

// initialize RingBuffer with default constructor function.
template<typename TP>
RingBuffer<TP>::RingBuffer() : _M_capacity(RING_BUFFER_INIT_SIZE) {
    _M_content = new TP[_M_capacity];
// set read/write position to 0;
    _M_read_pos = _M_write_pos = 0;
}

// destroy RingBuffer with de-constructor function.
template<typename TP>
RingBuffer<TP>::~RingBuffer() {
    delete[] _M_content;
}

// get the capacity of ringbuffer.
template<typename TP>
size_t RingBuffer<TP>::get_capacity() const {
    return _M_capacity;
}

// check whether the ringbuffer is empty
template<typename TP>
bool RingBuffer<TP>::is_empty() const {
// if the two pointers are the same, then the buffer is empty
    return _M_read_pos == _M_write_pos;
}

// read in one element from data
template<typename TP>
bool RingBuffer<TP>::read(TP &data) {
    // if the buffer is empty, then return false
    if (is_empty()) return false;
    // set the read data to data.
    data = _M_content[_M_read_pos % _M_capacity];
    // update the read pointer.
    _M_read_pos++;
    return true;
}


template<typename TP>
void RingBuffer<TP>::grow() {
    // tmp_1 is the new content, tmp_2 is the original one.
    TP *tmp_1, *tmp_2 = _M_content;
    if (_M_capacity < 1024) {
        tmp_1 = new TP[(int) (RING_BUFFER_GROW_FACTOR1 * _M_capacity)];
        //new content
    } else {
        tmp_1 = new TP[(int) (RING_BUFFER_GROW_FACTOR2 * _M_capacity)];
        //new content
    }
    // transfer value in content to tmp_1
    for (size_t i = 0; i < _M_capacity - 1; ++i) {
        read(tmp_1[i]);
    }
    //upd read pointer to 0
    _M_read_pos = 0;
    //upd write pointer ot capacity - 1
    _M_write_pos = _M_capacity - 1;
    // upd capacity
    _M_capacity *= (_M_capacity < 1024 ? RING_BUFFER_GROW_FACTOR1 : RING_BUFFER_GROW_FACTOR2);
    // upd content
    _M_content = tmp_1;
    // free the allocated memory
    delete[] tmp_2;

}


template<typename TP>
bool RingBuffer<TP>::write(const TP &data) {
    if (_M_write_pos + 1 - _M_read_pos == _M_capacity) {
        // if buffer will be full after writing, we upd capacity
        grow();
        // write the data and upd write pointer
        _M_content[_M_write_pos++] = data;
        return true;
    }
    // if buffer won't be full after insertion, just upd write pointer and write in data.
    _M_content[_M_write_pos % _M_capacity] = data;
    _M_write_pos++;
    return true;
}

// read multiple data via read function
template<typename TP>
bool RingBuffer<TP>::read_multi(size_t rdsize, std::vector<TP> &data) {
    // if there are no enough spaces to read, return false.
    if (_M_write_pos - _M_read_pos < rdsize) return false;
    // resize data
    if (rdsize > data.size()) {
        data.resize(rdsize, 0);
    }
    // read data via read function
    for (size_t i = 0; i < rdsize; ++i) {
        if (!read(data[i])) return false;
    }
    return true;
}

// write multiple data via write function
template<typename TP>
bool RingBuffer<TP>::write_multi(size_t wrtsize, const std::vector<TP> &data) {
    if (wrtsize > data.size()) return false;
    for (size_t i = 0; i < wrtsize; ++i) {
        // once failed in write a single data, we return false.
        if (!write(data[i])) return false;
    }
    return true;
}

// map function
template<typename TP>
bool RingBuffer<TP>::map(std::function<TP(TP)> &&func) {
    for (size_t i = _M_read_pos % _M_capacity; i < _M_write_pos % _M_capacity; ++i) {
        // apply func to every element.
        _M_content[i] = func(_M_content[i]);
    }
    return true;
}

// iterator class
template<typename TP>
class __detail::__iterator {
public:
    // define types
    // RB
    using RB = RingBuffer<TP>;
    // pointer
    using pointer = TP *;
    // reference
    using reference = TP &;
    // difference_type
    using difference_type = ptrdiff_t;
    // self_type
    using self_type = __detail::__iterator<TP>;

    // default constructor
    __iterator() = default;

    // copy constructor
    __iterator(const self_type &that) : ptr(that.ptr), buffer(that.buffer) {}

    // construct by two pointers
    __iterator(pointer that, RB *Buffer) : ptr(that), buffer(Buffer) {};

    // de-constructor
    ~__iterator() = default;

    // operator =
    self_type &operator=(self_type that) {
        this->ptr = that.ptr;
        return *this;
    }

    // deference *
    reference operator*() { return *ptr; }

    // operator ->
    pointer operator->() { return ptr; }

    // operator ==
    bool operator==(const self_type that) { return this->ptr == that.ptr; }

    // operator !=
    bool operator!=(const self_type that) { return this->ptr != that.ptr; }

    // prefix ++
    self_type &operator++() {
        ptr++;
        // if ptr exceed the boundary, we subtract the pointer with capacity
        if (ptr >= buffer->_M_content + buffer->_M_capacity) {
            ptr -= buffer->_M_capacity;
        }
        return *this;
    }

    // prefix --
    self_type &operator--() {
        ptr--;
        // if ptr exceed the boundary, we add the pointer with capacity
        if (ptr < buffer->_M_content) {
            ptr += buffer->_M_capacity;
        }
        return *this;
    }

    // postfix ++
    self_type operator++(int) {
        // origin is the original pointer
        self_type origin = *this;
        ptr++;
        // if ptr exceed the boundary, we subtract the pointer with capacity
        if (ptr >= buffer->_M_content + buffer->_M_capacity) {
            ptr -= buffer->_M_capacity;
        }
        return origin;
    }

    // postfix --
    self_type operator--(int) {
        // origin is the original pointer
        self_type origin = *this;
        ptr--;
        // if ptr exceed the boundary, we add the pointer with capacity
        if (ptr < buffer->_M_content) {
            ptr += buffer->_M_capacity;
        }
        return origin;
    }

    // operator +
    self_type operator+(difference_type num) {
        // tmp is the impermanent variable
        TP *tmp = ptr;
        tmp += num;
        // if tmp exceed the boundary, subtract it with capacity
        while (tmp >= buffer->_M_content + buffer->_M_capacity) {
            tmp -= buffer->_M_capacity;
        }
        return self_type(tmp, buffer);
    }

    // operator -
    self_type operator-(difference_type num) {
        // tmp is the impermanent variable
        TP *tmp = ptr;
        tmp -= num;
        // if tmp exceed the boundary, add it with capacity
        while (tmp < buffer->_M_content) {
            tmp += buffer->_M_capacity;
        }
        return self_type(tmp, buffer);
    }

    // operator +=
    self_type operator+=(difference_type num) {
        if (num < 0) {
            // if num is negative, then call -=
            return (*this) -= num;
        }
        ptr += num;
        // if ptr exceed the boundary, subtract the pointer with capacity
        while (ptr >= buffer->_M_content + buffer->_M_capacity) {
            ptr -= buffer->_M_capacity;
        }
        return *this;
    }

    // operator -=
    self_type operator-=(difference_type num) {
        if (num < 0) {
            // if num is negative, then call +=
            return (*this) += num;
        }
        ptr -= num;
        // if ptr exceed the boundary, add the pointer with capacity
        while (ptr < buffer->_M_content) {
            ptr += buffer->_M_capacity;
        }
        return *this;
    }


    // ptr is the pointer pointing to the element we are iterating
    TP *ptr;
    // buffer is the whole array holds the elements.
    RB *buffer;
};

// const_iterator class
template<typename TP>
class __detail::__const_iterator {
public:
    // define types
    // RB
    using RB = RingBuffer<TP>;
    // pointer
    using pointer = TP *;
    // reference
    using reference = TP &;
    // difference
    using difference_type = ptrdiff_t;
    // self_type
    using self_type = __detail::__const_iterator<TP>;

    // default constructor
    __const_iterator() = default;

    // copy constructor
    __const_iterator(const __iterator <TP> &that) : ptr(that.ptr), buffer(that.buffer) {}

    __const_iterator(self_type const &that) : ptr(that.ptr), buffer(that.buffer) {}

    // construct by two pointers
    __const_iterator(pointer const &that, const RB *Buffer) : ptr(that), buffer((RB *) Buffer) {}


// de-constructor
    ~__const_iterator() = default;

    // operator =
    self_type &operator=(self_type that) {
        this->ptr = that.ptr;
        this->buffer = that.buffer;
        return *this;
    }

    // deference *
    reference operator*() const { return *ptr; }

    // operator ->
    pointer operator->() { return ptr; }

    // operator ==
    bool operator==(const self_type that) { return this->ptr == that.ptr; }

    // operator !=
    bool operator!=(const self_type that) { return this->ptr != that.ptr; }

    // prefix ++
    self_type &operator++() {
        ptr++;
        // if ptr exceed the boundary, we subtract the pointer with capacity
        if (ptr >= buffer->_M_content + buffer->_M_capacity) {
            ptr -= buffer->_M_capacity;
        }
        return *this;
    }

    // prefix --
    self_type &operator--() {
        ptr--;
        // if ptr exceed the boundary, we add the pointer with capacity
        if (ptr < buffer->_M_content) {
            ptr += buffer->_M_capacity;
        }
        return *this;
    }

    // postfix ++
    self_type operator++(int) {
        // origin is the original pointer
        self_type origin = *this;
        ptr++;
        // if ptr exceed the boundary, we subtract the pointer with capacity
        if (ptr >= buffer->_M_content + buffer->_M_capacity) {
            ptr -= buffer->_M_capacity;
        }
        return origin;
    }

    // postfix --
    self_type operator--(int) {
        // origin is the original pointer
        self_type origin = *this;
        ptr--;
        // if ptr exceed the boundary, we add the pointer with capacity
        if (ptr < buffer->_M_content) {
            ptr += buffer->_M_capacity;
        }
        return origin;
    }

    // operator +
    self_type operator+(difference_type num) {
        // tmp is the impermanent variable
        TP *tmp = ptr;
        tmp += num;
        // if tmp exceed the boundary, subtract it with capacity
        while (tmp >= buffer->_M_content + buffer->_M_capacity) {
            tmp -= buffer->_M_capacity;
        }
        return self_type(tmp, buffer);
    }

    // operator -
    self_type operator-(difference_type num) {
        // tmp is the impermanent variable
        TP *tmp = ptr;
        tmp -= num;
        // if tmp exceed the boundary, add it with capacity
        while (tmp < buffer->_M_content) {
            tmp += buffer->_M_capacity;
        }
        return self_type(tmp, buffer);
    }

    // operator +=
    self_type operator+=(difference_type num) {
        if (num < 0) {
            // if num is negative, then call -=
            return (*this) -= -num;
        }
        ptr += num;
        // if ptr exceed the boundary, subtract the pointer with capacity
        while (ptr >= buffer->_M_content + buffer->_M_capacity) {
            ptr -= buffer->_M_capacity;
        }
        return *this;
    }

    // operator -=
    self_type operator-=(difference_type num) {
        if (num < 0) {
            // if num is negative, then call +=
            return (*this) += -num;
        }
        ptr -= num;
        // if ptr exceed the boundary, add the pointer with capacity
        while (ptr < buffer->_M_content) {
            ptr += buffer->_M_capacity;
        }
        return *this;
    }

    // ptr is the pointer pointing to the element we are iterating
    pointer ptr;
    // buffer is the whole array holds the elements.
    RB *buffer;
};

// begin
template<typename TP>
typename RingBuffer<TP>::iterator RingBuffer<TP>::begin() {
    // return the read_pos element iterator
    return iterator(&(_M_content[_M_read_pos % _M_capacity]), this);
}


// end
template<typename TP>
typename RingBuffer<TP>::iterator RingBuffer<TP>::end() {
    // return the write_pos element iterator
    return iterator(&(_M_content[_M_write_pos % _M_capacity]), this);
}


// cbegin
template<typename TP>
typename RingBuffer<TP>::const_iterator RingBuffer<TP>::cbegin() const {
    // return the read_pos element const iterator
    return const_iterator(&(_M_content[_M_read_pos % _M_capacity]), this);
}


// cend
template<typename TP>
typename RingBuffer<TP>::const_iterator RingBuffer<TP>::cend() const {
    // return the write_pos element const iterator
    return const_iterator(&(_M_content[_M_write_pos % _M_capacity]), this);
}










