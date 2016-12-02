#include <iostream>
#include <assert.h>
#include <array>
#include <cstdlib>
#include <sys/time.h>

template <typename _T, int _d>
struct BTreeNode
{
	static constexpr int d = _d;
	static constexpr int N = 2*d+1;
	using T = _T;
	using Pointer = BTreeNode*;

	BTreeNode(T t)
	{
		k = 1;
		data[0] = t;
		next[0] = next[1] = nullptr;
	}

	template <typename Iter>
	BTreeNode(Iter begin, Iter end)
	{
		k = 0;
		while(begin != end)
		{
			next[k] = nullptr;
			data[k++] = *begin++;
		}
		next[k] = nullptr;
	}

	int              k;
	std::array<T,N>     data;
	std::array<Pointer,N+1> next;
};

template <typename Node> using Pointer = typename Node::Pointer;
template <typename Node> using Data = typename Node::T;

template <typename Node>
Pointer<Node> split(Pointer<Node> head, Data<Node>& median)
{
	assert(head->k == N);

	Pointer<Node> new_right = new Node(head->data.begin() + Node::d + 1, head->data.end());
	head->k = Node::d;

	if(head->next[0] == nullptr)
	{
		new_right->next[0] = nullptr;
	}
	else
	{
		for(int i = 0; i <= Node::d; ++i)
		{
			new_right->next[i] = head->next[Node::d + i + 1];
		}
	}

	median = head->data[Node::d];
	return new_right;
}


template <typename Node>
bool insert(Pointer<Node> head, typename Node::T data, Data<Node>& up, Pointer<Node>& right)
{
	if(head->next[0] == nullptr)
	{
		const auto k = head->k;

		int i = 0;
		while(i < k && head->data[i] < data)
		{
			++i;
		}
		for(int j = k-1; j >= i; --j)
		{
			head->data[j+1] = head->data[j];
		}
		head->data[i] = data;
		head->k = k+1;
	}
	else
	{
		const auto k = head->k;

		int i = 0;
		while(i < k && head->data[i] < data)
			++i;

		Data<Node> tmp_up;
		Pointer<Node> tmp_right;
		if(insert<Node>(head->next[i], data, tmp_up, tmp_right))
		{
			for(int j = k-1; j >= i; --j)
			{
				head->data[j+1] = head->data[j];
				head->next[j+2] = head->next[j+1];
			}
			head->data[i] = tmp_up;
			head->next[i+1] = tmp_right;
			head->k = k+1;
		}
	}

	if(head->k == Node::N)
	{
		right = split<Node>(head, up);
		return true;
	}
	else
	{
		return false;
	}
}

template <typename Node>
Pointer<Node> insert(Pointer<Node> head, Data<Node> data)
{
	if(head == nullptr)
	{
		return new Node(data);
	}
	else
	{
		Pointer<Node> new_right;
		Data<Node> up;
		if(insert<Node>(head, data, up, new_right))
		{
			Pointer<Node> nn = new Node(up);
			nn->next[0] = head;
			nn->next[1] = new_right;
			return nn;
		}
		else
		{
			return head;
		}
	}
}

template <typename Node>
bool get(Pointer<Node> head, Data<Node> data)
{
	if(head->next[0] == nullptr)
	{
		for(int i = 0; i < head->k; ++i)
		{
			if(data == head->data[i])
			{
				return true;
			}
		}
		return false;
	}
	else
	{
		for(int i = 0; i < head->k; ++i)
		{
			if(data < head->data[i])
			{
				return get<Node>(head->next[i], data);
			}
			else if(data == head->data[i])
			{
				return true;
			}
		}
		return get<Node>(head->next[head->k], data);
	}
}

template <typename Node>
Pointer<Node> rebalance(Pointer<Node> head, int i)
{
	if(head->next[i]->k < Node::d)
	{
		if(i > 0 && head->next[i-1]->k > Node::d)
		{
			Pointer<Node> left  = head->next[i-1];
			Pointer<Node> right = head->next[i];

			if(right->next[0] != nullptr)
				right->next[right->k + 1] = right->next[right->k];
			for(int j = right->k; j > 0; --j)
			{
				std::cout << "  : " << right->data[j-1] << std::endl;
				right->data[j] = right->data[j-1];
				if(left->next[0] != nullptr)
					right->next[j] = right->next[j-1];
			}
			right->data[0] = head->data[i-1];
			if(left->next[0] != nullptr)
				right->next[0] = left->next[left->k];
			++(right->k);

			head->data[i-1] = left->data[left->k-1];

			--(left->k);

			std::cout << "Rotate Right" << std::endl;
		}
		else if(i < head->k && head->next[i+1]->k > Node::d)
		{
			Pointer<Node> left  = head->next[i];
			Pointer<Node> right = head->next[i+1];

			left->data[left->k] = head->data[i];
			++(left->k);
			if(left->next[0] != nullptr)
				left->next[left->k] = right->next[0];

			head->data[i] = right->data[0];
			for(int j = 0; j < right->k - 1; ++j)
			{
				right->data[j] = right->data[j+1];
			}
			--(right->k);

			std::cout << "Rotate Left" << std::endl;
		}
		else
		{
			if(i < head->k)
			{
				Pointer<Node> left = head->next[i];
				Pointer<Node> right = head->next[i+1];

				left->data[(left->k)++] = head->data[i];
				for(int j = 0; j < right->k; ++j)
				{
					left->data[left->k] = right->data[j];
					if(left->next[0] != nullptr)
						left->next[left->k] = right->next[j];
					++(left->k);
				}
				if(left->next[0] != nullptr)
					left->next[left->k] = right->next[right->k];

				delete right;

				--(head->k);
				for(int j = i; j < head->k; ++j)
				{
					head->data[j] = head->data[j+1];
					head->next[j+1] = head->next[j+2];
				}
				std::cout << "Collapse Right" << std::endl;
			}
			else
			{
				Pointer<Node> left = head->next[i-1];
				Pointer<Node> right = head->next[i];

				left->data[left->k] = head->data[i-1];
				++(left->k);
				for(int j = 0; j < right->k; ++j)
				{
					left->data[left->k] = right->data[j];
					if(left->next[0] != nullptr)
						left->next[left->k] = right->next[j];
					++(left->k);
				}
				if(left->next[0] != nullptr)
					left->next[left->k] = right->next[right->k];

				delete right;
				--(head->k);
				std::cout << "Collapse Left" << std::endl;
			}
		}
	}

//	std::cout << 
/*
	if(head->k == 0)
	{
		Pointer<Node> rval = head->next[0];
		std::cout << "Delete 3" << std::endl;
		delete head;
		std::cout << "Delete 3" << std::endl;
		return rval;
	}
	else
	{
		return head;
	}
*/
	return head;
}

template <typename Node>
void get_replacement(Pointer<Node> head, Data<Node>& key)
{
	if(head->next[0] == nullptr)
	{
		key = head->data[head->k-1];
		--(head->k);
	}
	else
	{
		get_replacement<Node>(head->next[head->k], key);
		rebalance<Node>(head, head->k);
	}
}

template <typename Node>
Pointer<Node> remove_H(Pointer<Node> head, Data<Node> data)
{
	if(head->next[0] == nullptr)
	{
		for(int i = 0; i < head->k; ++i)
		{
			if(data == head->data[i])
			{
				for(int j = i+1; j < head->k; ++j)
				{
					head->data[j-1] = head->data[j];
				}
				--(head->k);
				break;
			}
		}
		return head;
	}
	else
	{
		for(int i = 0; i < head->k; ++i)
		{
			if(data < head->data[i])
			{
				remove_H<Node>(head->next[i], data);
				return rebalance<Node>(head, i);
			}
			else if(data == head->data[i])
			{
				get_replacement<Node>(head->next[i], head->data[i]);
				return rebalance<Node>(head, i);
			}
		}
		remove_H<Node>(head->next[head->k], data);
		return rebalance<Node>(head, head->k);
	}
}

template <typename Node>
Pointer<Node> remove(Pointer<Node> head, Data<Node> data)
{
	head = remove_H<Node>(head, data);

	if(head->k == 0)
	{
		Pointer<Node> rval = head->next[0];
		std::cout << "Delete 3" << std::endl;
		delete head;
		std::cout << "Delete 3" << std::endl;
		return rval;
	}
	else
	{
		return head;
	}
}

template <typename T, int d>
std::ostream& operator<<(std::ostream& out, const BTreeNode<T,d>& head)
{
	out << "[ ";
	for(int i = 0; i < head.k; ++i)
	{
		if(head.next[0] != nullptr)
			out << *(head.next[i]) << " ";
		out << head.data[i] << " ";
	}
	if(head.next[0] != nullptr)
		out << *(head.next[head.k]);
	out << "]";
	return out;
}

template <typename T>
struct Interval
{
	Interval(T a, T b) : _a(a), _b(b) { assert(a <= b); }

	bool has(T x) { return _a <= x && x < _b; }

	T lower() const { return _a; }
	T upper()   const { return _b; }

private:
	T _a;
	T _b;
};

template <typename T>
bool operator<(const Interval<T>& x, const Interval<T>& y)
{
	return x.upper() <= y.lower();
}

template <typename T>
bool operator>(const Interval<T>& x, const Interval<T>& y)
{
	return x.lower() >= y.upper();
}

using BT = BTreeNode<int,1>;

int main(int argc, char *argv[])
{
	auto A = Interval<int>(0,2);
	auto B = Interval<int>(3,5);
	auto C = Interval<int>(2,3);

	std::cout << (A < B) << std::endl;
	std::cout << (A < C) << std::endl;
	std::cout << (C < B) << std::endl;
	std::cout << (B > A) << std::endl;

	std::cout << (A > B) << std::endl;
	std::cout << (B < A) << std::endl;

	for(auto i = A.lower(); i < A.upper(); ++i)
	{
		std::cout << i << " ";
	}
	std::cout << std::endl;

   timespec ts{0,0};
   clock_settime(CLOCK_THREAD_CPUTIME_ID, &ts);

	Pointer<BT> head = nullptr;
	for(int i = 0; i < 30000000; ++i)
	{
		auto add = std::rand();
		head = insert<BT>(head, add);
	}

   clock_gettime(CLOCK_THREAD_CPUTIME_ID, &ts); // Works on Linux

   std::cout << (ts.tv_sec * 10e9 + ts.tv_nsec)/10e9 << std::endl;

	int subtr[] = {49, 23, 58, 72, 44, 30, 78, 73, 87, 40, 29, 40, 42, 12, 92, 27, 9, 3, 7, 3, 9, 57, 33, 35, 16, 99, 78, 65, 60};

	for(int x : subtr)
	{
		head = remove<BT>(head, x);
//		std::cout << x << " - " << *head << std::endl;
	}

//	std::cout << *head << std::endl;

/*
	Data<BT> x;
	std::cout << " 1: ";
	print<BT>(split<BT>(head, x));
	std::cout << x << std::endl;
	std::cout << " 2: ";
	print<BT>(head);
	*/
}
