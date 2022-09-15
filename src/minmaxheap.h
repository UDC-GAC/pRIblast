/*
 Copyright (C) 2018 Kilian Gebhardt

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#ifndef MINMAXHEAP_H
#define MINMAXHEAP_H

#include <vector>

#include <stdint.h>

// see https://stackoverflow.com/a/24748637
inline int uint64_log2(uint64_t n)
{
  #define S(k) if (n >= (UINT64_C(1) << k)) { i += k; n >>= k; }

  int i = -(n == 0); S(32); S(16); S(8); S(4); S(2); S(1); return i;

  #undef S
}

namespace minmaxheap {

	template <typename T>
	class MinMaxHeap {
		std::vector<T> heap;
		private:
			void trickledown(const size_t i);
			inline void trickledownmin(size_t i);
			inline void trickledownmax(size_t i);
			void bubbleup(const size_t i);
			inline void bubbleupmin(size_t i);
			inline void bubbleupmax(size_t i);
			void swap(const size_t i, const size_t j);
			inline short level(const size_t) const;

		public:
			MinMaxHeap(size_t reserve=0);
			void reserve(size_t n);
			inline size_t size() const;
			inline void clear();
			void insert(T key);
			T peekmin() const;
			T peekmax() const;
			T popmin();
			T popmax();
			std::vector<T> getheap() const {
				return heap;
			}
	};

	template<typename T>
	MinMaxHeap<T>::MinMaxHeap(size_t reserve) {
		heap.reserve(reserve);
	}

	template<typename T>
	void MinMaxHeap<T>::reserve(size_t n) {
		heap.reserve(n);
	}

	template<typename T>
	inline size_t MinMaxHeap<T>::size() const {
		return heap.size();
	}

	template<typename T>
	inline void MinMaxHeap<T>::clear() {
		heap.clear();
	}

	template<typename T>
	void MinMaxHeap<T>::insert(T key) {
		heap.push_back(key);
		bubbleup(heap.size() - 1);
	}

	template<typename T>
	T MinMaxHeap<T>::peekmax() const {
		assert(heap.size() > 0);
		if (heap.size() == 1)
			return heap[0];
		else if (heap.size() == 2)
			return heap[1];
		else
			return heap[1] > heap[2] ? heap[1] : heap[2];
	}

	template<typename T>
	T MinMaxHeap<T>::popmax() {
		T e;
		if (heap.size() == 1) {
			e = heap[0];
			heap.pop_back();
		}
		else if (heap.size() == 2) {
			e = heap[1];
			heap.pop_back();
		}
		else {
			const size_t i {heap[1] > heap[2] ? (size_t) 1 : (size_t) 2};
			e = heap[i];
			heap[i] = heap.back();
			heap.pop_back();
			trickledown(i);
		}
		return e;
	}

	template<typename T>
	T MinMaxHeap<T>::peekmin() const {
		return heap.front();
	}

	template<typename T>
	T MinMaxHeap<T>::popmin() {
		const T e {heap[0]};
		heap[0] = heap.back();
		heap.pop_back();
		trickledown(0);
		return e;
	}

	template<typename T>
	void MinMaxHeap<T>::trickledown(size_t i) {
		if (level(i) % 2 == 0)  // min level
			trickledownmin(i);
		else
			trickledownmax(i);
	}

	template<typename T>
	inline void MinMaxHeap<T>::trickledownmin(size_t i) {
		bool child {true};
		while (heap.size() > 2  * i + 1) {
			size_t m = 2 * i + 1;
			if (m+1 < heap.size() and heap[m+1] < heap[m])
				 ++m;
			const size_t bound {4 * i + 7 < heap.size()
				? 4 * i + 7 : heap.size()};
			for (size_t j = 4 * i + 3; j < bound; ++j)
				if (heap[j] < heap[m]) {
					m = j;
					child = false;
				}
			if (child) {
				if (heap[m] < heap[i])
					swap(m, i);
				break;
			} else {
				if (heap[m] < heap[i]) {
					swap(i, m);
					if (heap[m] > heap[(m-1) / 2])
						swap(m, (m-1) / 2);

//					trickledownmin(m);
					i = m;
					child = true;
				} else
					break;
			}
		}
	}

	template<typename T>
	inline void MinMaxHeap<T>::trickledownmax(size_t i) {
		bool child {true};
		while (heap.size() > 2  * i + 1) {
			size_t m = 2 * i + 1;
			if (m + 1 < heap.size() and heap[m+1] > heap[m])
				++m;
			const size_t bound {4 * i + 7 < heap.size()
				? 4 * i + 7 : heap.size()};
			for (size_t j = 4 * i + 3; j < bound; ++j)
				if (heap[j] > heap[m]) {
					m = j;
					child = false;
				}
			if (child) {
				if (heap[m] > heap[i])
					swap(m, i);
				break;
			}
			else {
				if (heap[m] > heap[i]) {
					swap(i,m);
					if (heap[m] < heap[(m-1) / 2])
						swap(m, (m-1) / 2);

//					trickledownmax(m);
					i = m;
					child = true;
				} else
					break;
			}
		}
	}

	template<typename T>
	void MinMaxHeap<T>::bubbleup(const size_t i) {
		if (level(i) % 2 == 0) { // min level
			if (i > 0 and heap[i] > heap[(i-1) / 2]) {
				swap(i, (i-1)/2);
				bubbleupmax((i-1) / 2);
			} else {
				bubbleupmin(i);
			}
		} else {
			if (i > 0 and heap[i] < heap[(i-1) / 2]) {
				swap(i, (i-1) / 2);
				bubbleupmin((i-1) / 2);
			} else {
				bubbleupmax(i);
			}
		}
	}

	template<typename T>
	inline void MinMaxHeap<T>::bubbleupmin(size_t i) {
		while (i > 2)
			if (heap[i] < heap[(i-3) / 4]) {
				swap(i, (i-3)/4);
				i = (i-3) / 4;
			} else
				break;
	}

	template<typename T>
	inline void MinMaxHeap<T>::bubbleupmax(size_t i) {
		while (i > 2)
			if (heap[i] > heap[(i-3) / 4]) {
				swap(i, (i-3)/ 4);
				i = (i-3) / 4;
			} else
				break;
	}

	template<typename T>
	inline short MinMaxHeap<T>::level(size_t i) const {
		return (uint64_log2(i + 1));
	}

	template<typename T>
	inline void MinMaxHeap<T>::swap(size_t i, size_t j) {
		const T tmp {heap[i]};
		heap[i] = heap[j];
		heap[j] = tmp;
	}
}

#endif
