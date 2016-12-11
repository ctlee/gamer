#include <iostream>
#include <assert.h>
#include <array>
#include <vector>
#include <cstdlib>
#include <limits>

template <typename T>
struct Interval
{
	Interval() : _a(0), _b(0) {}
	Interval(T a) : _a(a), _b(a+1) {}
	Interval(T a, T b) : _a(a), _b(b) { assert(a <= b); }
	Interval(const Interval<T>& rhs) : _a(rhs._a), _b(rhs._b) {}

	Interval& operator=(const Interval& rhs)
	{
		_a = rhs._a;
		_b = rhs._b;
		return *this;
	}

	bool has(T x) { return _a <= x && x < _b; }

	T  lower() const { return _a; }
	T  upper() const { return _b; }

	T& lower()       { return _a; }
	T& upper()       { return _b; }

	size_t size() { return _b - _a; }

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

template <typename T>
bool operator<(T x, const Interval<T>& y)
{
	return x < y.lower();
}

template <typename T>
bool operator>(const Interval<T>& x, T y)
{
	return x.lower() > y;
}

template <typename T>
bool operator<(const Interval<T>& x, T y)
{
	return x.upper() <= y;
}

template <typename T>
bool operator>(T x, const Interval<T>& y)
{
	return x >= y.upper();
}

template <typename T>
bool operator==(const Interval<T>& x, const Interval<T>& y)
{
	return (x.lower() == y.lower()) && (x.upper() && y.upper());
}

template <typename T>
std::ostream& operator<<(std::ostream& out, const Interval<T>& x)
{
	out << "[" << x.lower() << "~" << x.upper() << ")";
	return out;
}

template <typename T>
int merge(Interval<T>& A, T x)
{
	if(x + 1 < A.lower())
		return 0;
	else if(x + 1 == A.lower())
	{
		A.lower() = x;
		return 1;
	}
	else if(A.lower() <= x && x < A.upper())
		return 2;
	else if(A.upper() == x)
	{
		A.upper() = x + 1;
		return 3;
	}
	else if(A.upper() < x)
		return 4;
	else
		return 5;
}


template <typename _T, int _d>
struct BTreeNode
{
	static constexpr int d = _d;
	static constexpr int N = 2*d+1;
	using Scalar = _T;
	using Data = Interval<Scalar>;
	using Pointer = BTreeNode*;

	BTreeNode() {}
	BTreeNode(const Data& t)
		: k(1), next{nullptr,nullptr}
	{
		data[0] = t;
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

	int                     k;
	std::array<Data,N>      data;
	std::array<Pointer,N+1> next;
};

template <typename Node> using Pointer = typename Node::Pointer;
template <typename Node> using Data = typename Node::Data;
template <typename Node> using Scalar = typename Node::Scalar;



template <typename Node>
void rebalance(Pointer<Node> head, int i)
{
	Pointer<Node> curr = head->next[i];

	if(curr->k == Node::N)
	{
		Pointer<Node> left = curr;
		Pointer<Node> right = new Node(curr->data.begin() + Node::d + 1, curr->data.end());
		curr->k = Node::d;

		if(curr->next[0] == nullptr)
		{
			right->next[0] = nullptr;
		}
		else
		{
			for(int i = 0; i <= Node::d; ++i)
			{
				right->next[i] = curr->next[Node::d + i + 1];
			}
		}

		Data<Node> up = curr->data[Node::d];

		for(int j = head->k - 1; j >= i; --j)
		{
			head->data[j+1] = head->data[j];
			head->next[j+2] = head->next[j+1];
		}
		head->data[i] = up;
		head->next[i+1] = right;
		++(head->k);
	}
	else if(curr->k < Node::d)
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
}

template <typename Node>
void insert_H(Pointer<Node> head, const Data<Node>& data)
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

		insert_H<Node>(head->next[i], data);
		rebalance<Node>(head, i);
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
		insert_H<Node>(head, data);
		if(head->k == Node::N)
		{
			Pointer<Node> nn = new Node();
			nn->k = 0;
			nn->next[0] = head;
			rebalance<Node>(nn, 0);
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
void remove_H(Pointer<Node> head, Data<Node> data)
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
	}
	else
	{
		for(int i = 0; i < head->k; ++i)
		{
			if(data < head->data[i])
			{
				remove_H<Node>(head->next[i], data);
				rebalance<Node>(head, i);
				return;
			}
			else if(data == head->data[i])
			{
				get_replacement<Node>(head->next[i], head->data[i]);
				rebalance<Node>(head, i);
				return;
			}
		}
		remove_H<Node>(head->next[head->k], data);
		rebalance<Node>(head, head->k);
	}
}

template <typename Node>
Pointer<Node> remove(Pointer<Node> head, Data<Node> data)
{
	remove_H<Node>(head, data);

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
void fill_left(Pointer<Node> head, Data<Node>& x)
{
	if(head->next[0] == nullptr)
	{
		Data<Node>& left = head->data[head->k-1];
		if(left.upper() == x.lower())
		{
			x.lower() = left.lower();
			--(head->k);
		}
	}
	else
	{
		fill_left<Node>(head->next[head->k], x);
		rebalance<Node>(head, head->k);
	}
}


template <typename Node>
void fill_right(Pointer<Node> head, Data<Node>& x)
{
	if(head->next[0] == nullptr)
	{
		Data<Node>& right = head->data[0];
		if(right.lower() == x.upper())
		{
			x.upper() = right.upper();
			--(head->k);
			for(int i = 0; i < head->k; ++i)
			{
				head->data[i] = head->data[i+1];
			}
		}
	}
	else
	{
		fill_right<Node>(head->next[0], x);
		rebalance<Node>(head, 0);
	}
}


template <typename Node>
void insert_scalar_H(Pointer<Node> head, Scalar<Node> data)
{
	if(head->next[0] == nullptr)
	{
		const auto k = head->k;

		int i;
		for(i = 0; i < k; ++i)
		{
			Data<Node>& A = head->data[i];
			Scalar<Node> x = data;

			if(x + 1 < A.lower())
			{
				for(int j = k-1; j >= i; --j)
				{
					head->data[j+1] = head->data[j];
				}
				head->data[i] = data;
				++(head->k);
				return;
			}
			else if(x + 1 == A.lower())
			{
				A.lower() = x;
				return;
			}
			else if(A.lower() <= x && x < A.upper())
			{
				return;
			}
			else if(A.upper() == x)
			{
				if(i + 1 < k)
				{
					Data<Node>& B = head->data[i+1];
					if(x + 1 == B.lower())
					{
						A.upper() = B.upper();
						for(int j = i+1; j < k-1; ++j)
						{
							head->data[j] = head->data[j+1];
						}
						--(head->k);
					}
					else
					{
						A.upper() = x + 1;
					}
				}
				else
				{
					A.upper() = x + 1;
				}
				return;
			}
		}
		head->data[i] = data;
		++(head->k);
	}
	else
	{
		const auto k = head->k;

		int i;
		for(i = 0; i < k; ++i)
		{
			Data<Node>& A = head->data[i];
			Scalar<Node> x = data;

			if(x + 1 < A.lower())
			{
				insert_scalar_H<Node>(head->next[i], data);
				rebalance<Node>(head, i);
				return;
			}
			else if(x + 1 == A.lower())
			{
				A.lower() = x;
				fill_left<Node>(head->next[i], A);
				rebalance<Node>(head, i);
				return;
			}
			else if(A.lower() <= x && x < A.upper())
			{
				return;
			}
			else if(A.upper() == x)
			{
				A.upper() = x + 1;
				fill_right<Node>(head->next[i+1], A);
				rebalance<Node>(head, i+1);
				return;
			}
		}
		insert_scalar_H<Node>(head->next[i], data);
		rebalance<Node>(head, i);
	}
}

template <typename Node>
Pointer<Node> insert_scalar(Pointer<Node> head, Scalar<Node> data)
{
	if(head == nullptr)
	{
		return new Node(data);
	}
	else
	{
		insert_scalar_H<Node>(head, data);
		if(head->k == Node::N)
		{
			Pointer<Node> nn = new Node();
			nn->k = 0;
			nn->next[0] = head;
			rebalance<Node>(nn, 0);
			return nn;
		}
		else if(head->k == 0)
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
}

template <typename Node>
void insert_left(Pointer<Node> head, const Data<Node>& x)
{
	if(head->next[0] == nullptr)
	{
		head->data[head->k] = x;
		++(head->k);
	}
	else
	{
		insert_left<Node>(head->next[head->k], x);
		rebalance<Node>(head, head->k);
	}
}


template <typename Node>
bool remove_scalar_H(Pointer<Node> head, Scalar<Node> x)
{
	if(head->next[0] == nullptr)
	{
		const auto k = head->k;

		int i;
		for(i = 0; i < k; ++i)
		{
			Data<Node>& A = head->data[i];

			if(x < A.lower())
			{
//				std::cout << "if(x < A.lower())" << std::endl;
				return false;
			}
			else if(x == A.lower())
			{
//				std::cout << "if(x == A.lower())" << std::endl;
				if(x + 1 == A.upper())
				{
//					std::cout << "if(x + 1 == A.upper())" << std::endl;
					--(head->k);
					for(int j = i; j < head->k; ++j)
					{
						head->data[j] = head->data[j+1];
					}
					return true;
				}
				A.lower() = x + 1;
				return true;
			}
			else if(/*A.lower() < x &&*/ x + 1 < A.upper())
			{
//				std::cout << "x + 1 < A.upper()" << std::endl;
				for(int j = head->k; j > i; --j)
				{
					head->data[j] = head->data[j-1];
				}
				++(head->k);
				A.upper() = x;
				head->data[i+1].lower() = x + 1;
				return true;
			}
			else if(x + 1 == A.upper())
			{
//				std::cout << "x + 1 < A.upper()" << std::endl;
				A.upper() = x;
				return true;
			}
		}
		return false;
	}
	else
	{
		const auto k = head->k;

		int i;
		for(i = 0; i < k; ++i)
		{
			Data<Node>& A = head->data[i];

			if(x < A.lower())
			{
				bool rval = remove_scalar_H<Node>(head->next[i], x);
				rebalance<Node>(head, i);
				return rval;
			}
			else if(x == A.lower())
			{
				if(x + 1 == A.upper())
				{
					get_replacement<Node>(head->next[i], A);
					rebalance<Node>(head, i);
					return true;
				}
				A.lower() = x + 1;
				return true;
			}
			else if(/*A.lower() < x &&*/ x + 1 < A.upper())
			{
				Data<Node> B(A.lower(), x);
				A.lower() = x + 1;
				insert_left<Node>(head->next[i], B);
				rebalance<Node>(head, i);
				return true;
			}
			else if(x + 1 == A.upper())
			{
				A.upper() = x;
				return true;
			}
		}
		bool rval = remove_scalar_H<Node>(head->next[i], x);
		rebalance<Node>(head, i);
		return rval;
	}
}

template <typename Node>
bool remove_scalar(Pointer<Node>& head, Scalar<Node> data)
{
	if(head == nullptr)
	{
		return false;
	}

	bool rval = remove_scalar_H<Node>(head, data);

	if(head->k == Node::N)
	{
		Pointer<Node> nn = new Node();
		nn->k = 0;
		nn->next[0] = head;
		rebalance<Node>(nn, 0);
		head = nn;
	}
	else if(head->k == 0)
	{
		Pointer<Node> tmp = head;
		head = head->next[0];
		delete tmp;
	}

	return rval;
}

template <typename Node>
Scalar<Node> pop_scalar(Pointer<Node>& head)
{
	if(head)
	{
		Scalar<Node> x = head->data[0].lower();
		remove_scalar<Node>(head, x);
		return x;
	}
}

template <typename Node>
void destruct(Pointer<Node> head)
{
	if(head == nullptr)
	{
		return;
	}
	else
	{
		if(head->next[0] != nullptr)
		{
			for(int i = 0; i < head->k; ++i)
			{
				destruct<Node>(head->next[i]);
			}
		}
		delete head;
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

template <typename T, int d>
std::ostream& operator<<(std::ostream& out, const BTreeNode<T,d>* head)
{
	if(head == nullptr)
	{
		out << "[nil]";
	}
	else
	{
		out << "[ ";
		for(int i = 0; i < head->k; ++i)
		{
			if(head->next[0] != nullptr)
				out << head->next[i] << " ";
			out << head->data[i] << " ";
		}
		if(head->next[0] != nullptr)
			out << head->next[head->k];
		out << "]";
	}
	return out;
}

template <typename _T, int _d = 16>
class index_tracker
{
public:
	using Node = BTreeNode<_T, _d>;
	using T = _T;
	constexpr static int d = _d;

	index_tracker() : head(new Node(Interval<T>(0,std::numeric_limits<T>::max()))) {}
	~index_tracker()
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



	friend std::ostream& operator<<(std::ostream& out, const index_tracker& x)
	{
		out << x.head;
		return out;
	}

private:
	Pointer<Node> head;
};
