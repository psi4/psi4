/*
 * MIT License
 *
 * Copyright (c) 2018 Tobias Widlund
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#pragma once
#include <type_traits>
#include <cstddef>
#include <utility>

namespace en {
template <typename container_type>
struct enumerate_wrapper {
    using iterator_type = std::conditional_t<std::is_const_v<container_type>, typename container_type::const_iterator,
                                             typename container_type::iterator>;
    using pointer_type = std::conditional_t<std::is_const_v<container_type>, typename container_type::const_pointer,
                                            typename container_type::pointer>;
    using reference_type = std::conditional_t<std::is_const_v<container_type>, typename container_type::const_reference,
                                              typename container_type::reference>;

    constexpr enumerate_wrapper(container_type& c) : container(c) {}

    struct enumerate_wrapper_iter {
        size_t index;
        iterator_type value;

        constexpr bool operator!=(const iterator_type& other) const { return value != other; }
        constexpr enumerate_wrapper_iter& operator++() {
            ++index;
            ++value;
            return *this;
        }

        constexpr std::pair<size_t, reference_type> operator*() {
            return std::pair<size_t, reference_type>{index, *value};
        }
    };

    constexpr enumerate_wrapper_iter begin() { return {0, std::begin(container)}; }

    constexpr iterator_type end() { return std::end(container); }
    container_type& container;
};

template <typename container_type>
constexpr auto enumerate(container_type& c) {
    return enumerate_wrapper(c);
}

template <typename container_type>
constexpr auto const_enumerate(const container_type& c) {
    return enumerate_wrapper(c);
}
}  // namespace en
