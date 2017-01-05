#include <iostream>
#include <assert.h>
#include <array>
#include <vector>
#include <set>
#include <cstdlib>
#include <sys/time.h>
#include <stack>

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

//			std::cout << "Rotate Right" << std::endl;
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
				if(right->next[0] != nullptr)
					right->next[j] = right->next[j+1];
			}
			--(right->k);
			if(right->next[0] != nullptr)
				right->next[right->k] = right->next[right->k + 1];

//			std::cout << "Rotate Left" << std::endl;
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
		delete head;
		return rval;
	}
	else
	{
		return head;
	}
}


template <typename Node>
Data<Node> check_order(Pointer<Node> head, Data<Node> curr)
{
	if(head != nullptr)
	{		
		if(head->next[0] == nullptr)
		{
			for(int i = 0; i < head->k; ++i)
			{
				if(curr > head->data[i])
				{
					std::cout << "ORDER WRONG!!!   --   " << curr << " > " << head->data[i] << std::endl;
					exit(1);
				}
				curr = head->data[i];
			}
		}
		else
		{
			for(int i = 0; i < head->k; ++i)
			{
				curr = check_order<Node>(head->next[i], curr);
				if(curr > head->data[i])
				{
					std::cout << "ORDER WRONG!!!   --   " << curr << " > " << head->data[i] << std::endl;
					exit(1);
				}
				curr = head->data[i];
			}
			curr = check_order<Node>(head->next[head->k], curr);
		}
	}
	return curr;
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

template <typename Node>
struct btree_iterator : public std::iterator<std::bidirectional_iterator_tag, Data<Node>> {
public:
	using super = std::iterator<std::bidirectional_iterator_tag, Data<Node>>;
	btree_iterator() {}
	btree_iterator(const Pointer<Node> head)
	{
		push(head);
	}

	btree_iterator& operator++()
	{
		Pointer<Node> head = curr.top().first;
		int &i = curr.top().second;

		++i;
		if(head->next[0] != nullptr)
		{
			push(head->next[i]);
		}
		else if(i == head->k)
		{
			pop();
		}
/*
		if(!curr.empty() && curr.top().first->k == curr.top().second)
		{
			std::cout << "# " << i << " #" << std::endl;
		}
*/
		return *this;
	}

	btree_iterator operator++(int) { auto tmp = *this; ++(*this); return tmp; }

	btree_iterator& operator--()
	{
		return *this;
	}

	btree_iterator operator--(int) { auto tmp = *this; --(*this); return tmp; }

	bool operator==(const btree_iterator& j) const
	{
		return curr == j.curr;
	}

	bool operator!=(const btree_iterator& j) const { return !(*this == j); }

	Data<Node> operator*()
	{
		return curr.top().first->data[curr.top().second];
	}

    typename super::pointer operator->() { return curr.top().first->data[curr.top().second]; }

protected:
	void push(Pointer<Node> p)
	{
		while(p != nullptr)
		{
			curr.push(std::make_pair(p,0));
			p = p->next[0];
		}
	}

	void pop()
	{
		while(!curr.empty() && curr.top().first->k == curr.top().second)
		{
			curr.pop();
		}
	}

	std::stack<std::pair<Pointer<Node>,int>> curr;
};

template <typename Node>
bool check_order(Pointer<Node> head)
{
	Data<Node> last = std::numeric_limits<Data<Node>>::min();

	btree_iterator<Node> curr(head);
	btree_iterator<Node> end;

	while(curr != end)
	{
		if(*curr < last)
		{
			return false;
		}
		last = *curr;
		++curr;
	}
	return true;
}

template <typename T, int d>
std::ostream& operator<<(std::ostream& out, BTreeNode<T,d>* head)
{
	if(head != nullptr)
	{
		btree_iterator<BTreeNode<T,d>> curr(head);
		btree_iterator<BTreeNode<T,d>> end;

		while(curr != end)
		{
			out << *curr << " ";
			++curr;
		}
	}
	return out;
}
/*

template <typename _T, int _d = 16>
class btree_set
{
public:
	using Node = BTreeNode<_T, _d>;
	using T = _T;
	constexpr static int d = _d;

	btree_set() : head(new Node(Interval<T>(0,std::numeric_limits<T>::max()))) {}
	~btree_set()
	{
		destruct<Node>(head);
	}

	void insert(T x)
	{
		head = insert_scalar<Node>(head, x);
	}

	Scalar<Node> pop()
	{
		auto x = pop_scalar<Node>(head);
    	return x;
	}

	void remove(Scalar<Node> x)
	{
		remove_scalar<Node>(head, x);
	}

	bool empty() const
	{
		return head == nullptr;
	}

	friend std::ostream& operator<<(std::ostream& out, const btree_set& x)
	{
		out << x.head;
		return out;
	}

private:
	Pointer<Node> head;
};
*/

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
	std::set<int> benchmark;
	std::vector<int> added_list;
	constexpr size_t N = 100;
	for(int i = 0; i < N; ++i)
	{
		int add = std::rand() % 1000000;
		head = insert<BT>(head, i);
//		benchmark.insert(add);
//		check_order<BT>(head, -1);
		added_list.push_back(add);
		std::cout << head << std::endl;
//		std::cout << check_order<BT>(head) << std::endl;
	}

	clock_gettime(CLOCK_THREAD_CPUTIME_ID, &ts); // Works on Linux

	std::cout << (ts.tv_sec * 10e9 + ts.tv_nsec)/10e9 << std::endl;

	for(int x : added_list)
	{
//		check_order<BT>(head, -1);
		head = remove<BT>(head, x);
//		std::cout << x << " - " << *head << std::endl;
//		check_order<BT>(head, -1);
//		benchmark.erase(x);
	}

	if(head)
	{
			std::cout << head << std::endl;
	}
	else
	{
		std::cout << "Empty!" << std::endl;
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
