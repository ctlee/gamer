
template <class Data>
struct Node {};

template <class Data>
struct Root : Node<Data> {
	Node<Data>* gt;
	Node<Data>* lt;
};

template <class Data>
struct Internal : Node<Data> {
	Node<Data>* parent;
	Node<Data>* gt;
	Node<Data>* lt;
};

template <class Data>
struct Leaf : Node<Data> {
	Node<Data>* parent;
	Data data;
};

template <class Data>
struct kdtree {
	template <class Iterator>
	kdtree(Iterator begin, Iterator end)
	{
		std::size_t size = end - begin;
	}
};
