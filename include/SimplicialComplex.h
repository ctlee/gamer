#pragma once

#include <map>
#include <set>
#include <iterator>
#include <array>
#include <vector>
#include <iostream>
#include "util.h"

//using Simplex = unsigned int;

namespace detail {
	template <class T> using map = std::map<size_t,T>;
	/*
	 * asc_Node must be defined outside of simplicial_complex because c++ does not allow
	 * internal templates to be partially specialized. I admit that I do not understand
	 * the reason for this, but this work around seems to work. However, the muddying of
	 * the template parameters is unfortunate.
	 */
	template <class KeyType, size_t k, size_t N, typename DataTypes, class> struct asc_Node;

	/**
	 * @brief      A base class for a Simplicial Complex Node
	 */
	struct asc_NodeBase {
		/**
		 * @brief      Node constructor
		 *
		 * @param[in]  id    The identifier
		 */
		asc_NodeBase(int id) : _node(id) {}
		
		/**
		 * @brief      Destroys the object
		 */
		virtual ~asc_NodeBase() {};

		size_t _node;	/**< Internal Node ID*/
	};

	/**
	 * @brief      Base class for Node with some data
	 *
	 * @tparam     DataType  type of the data to be contained
	 */
	template <class DataType>
	struct asc_NodeData {
		DataType _data; /**< stored data with type DataType */
	};

	/**
	 * @brief      Base class for Node with DataType of void
	 */
	template <>
	struct asc_NodeData<void>
	{};

	/**
	 * @brief      Base class for Node with edge data of type DataType
	 *
	 * @tparam     KeyType   type for indexing Nodes
	 * @tparam     DataType  type of edge_data stored
	 */
	template <class KeyType, class DataType>
	struct asc_EdgeData
	{
		std::map<KeyType, DataType> _edge_data;
	};	

	/**
	 * @brief      Base class for Node with no edge data
	 *
	 * @tparam     KeyType  type for indexing Nodes
	 */
	template <class KeyType>
	struct asc_EdgeData<KeyType,void>
	{};

	/**
	 * @brief      Base class for Node with parent Nodes
	 *
	 * @tparam     KeyType        type for indexing Nodes
	 * @tparam     k              the depth of the Node
	 * @tparam     N              the maximum depth of Nodes in the simplicial complex
	 * @tparam     NodeDataTypes  a util::type_holder array of Node types
	 * @tparam     EdgeDataTypes  a util::type_holder array of Edge types
	 */
	template <class KeyType, size_t k, size_t N, class NodeDataTypes, class EdgeDataTypes>
	struct asc_NodeDown : public asc_EdgeData<KeyType,typename util::type_get<k-1,EdgeDataTypes>::type>
	{
		using DownNodeT = asc_Node<KeyType,k-1,N,NodeDataTypes,EdgeDataTypes>;
		std::map<KeyType, DownNodeT*> _down;	/**< @brief Map of pointers to parents */
	};

	/**
	 * @brief      Base class for Node with children Nodes
	 *
	 * @tparam     KeyType        type for indexing Nodes
	 * @tparam     k              the depth of the Node
	 * @tparam     N              the maximum depth of Nodes in the simplicial complex
	 * @tparam     NodeDataTypes  a util::type_holder array of Node types
	 * @tparam     EdgeDataTypes  a util::type_holder array of Edge types
	 */
	template <class KeyType, size_t k, size_t N, class NodeDataTypes, class EdgeDataTypes>
	struct asc_NodeUp {
		using UpNodeT = asc_Node<KeyType,k+1,N,NodeDataTypes,EdgeDataTypes>;
		std::map<KeyType, UpNodeT*> _up;	/**< @brief Map of pointers to children */
	};

	/**
	 * @brief      Node with both parents and children 
	 *
	 * @tparam     KeyType        index typename
	 * @tparam     k              level of the node
	 * @tparam     N              max level of the simplicial complex
	 * @tparam     NodeDataTypes  util::type_holder of node types
	 * @tparam     EdgeDataTypes  util::type_holder of edge types
	 */
	template <class KeyType, size_t k, size_t N, class NodeDataTypes, class EdgeDataTypes>
	struct asc_Node : public asc_NodeBase,
					  public asc_NodeData<typename util::type_get<k,NodeDataTypes>::type>,
					  public asc_NodeDown<KeyType, k, N, NodeDataTypes, EdgeDataTypes>,
					  public asc_NodeUp<KeyType, k, N, NodeDataTypes, EdgeDataTypes>
	{
		asc_Node(int id) : asc_NodeBase(id) {}
        
        friend std::ostream& operator<<(std::ostream& output, const asc_Node& node){
            output  << "Node(level=" << k << ", " << "id=" << node._node;
            if(node._down.size() > 0)
                for(auto it=node._down.cbegin(); it!=node._down.cend(); ++it)
                    output  << ", NodeDownID={'"
                            << it->first << "', "
                            << it->second->_node << "}";
            if(node._up.size() > 0)
               for(auto it=node._up.cbegin(); it!=node._up.cend(); ++it)
                    output  << ", NodeUpID={'"
                            << it->first << "', "
                            << it->second->_node << "}";
            output  << ")";
            return output;
        }
	};

	/**
	 * @brief      Root node with only children
	 *
	 * @tparam     KeyType        index typename
	 * @tparam     N              max level of the simplicial complex
	 * @tparam     NodeDataTypes  util::type_holder of node types
	 * @tparam     EdgeDataTypes  util::type_holder of edge types
	 */
	template <class KeyType, size_t N, class NodeDataTypes, class EdgeDataTypes>
	struct asc_Node<KeyType,0,N,NodeDataTypes,EdgeDataTypes> : public asc_NodeBase,
											 public asc_NodeData<typename util::type_get<0,NodeDataTypes>::type>,
											 public asc_NodeUp<KeyType, 0, N, NodeDataTypes, EdgeDataTypes>
	{
		asc_Node(int id) : asc_NodeBase(id) {}
       
        friend std::ostream& operator<<(std::ostream& output, const asc_Node& node){
            output  << "Node(level=" << 0
                    << ", id=" << node._node;
            if(node._up.size() > 0)
                for(auto it=node._up.cbegin(); it!=node._up.cend(); ++it)
                    output  << ", NodeUpID={'"
                            << it->first << "', "
                            << it->second->_node << "}";
            output << ")";
            return output;
        }
	};

	/**
	 * @brief      Top level node with only parents
	 *
	 * @tparam     KeyType        index typename
	 * @tparam     N              max level of the simplicial complex
	 * @tparam     NodeDataTypes  util::type_holder of node types
	 * @tparam     EdgeDataTypes  util::type_holder of edge types
	 */
	template <class KeyType, size_t N, class NodeDataTypes, class EdgeDataTypes>
	struct asc_Node<KeyType,N,N,NodeDataTypes,EdgeDataTypes> : public asc_NodeBase,
											 public asc_NodeData<typename util::type_get<N,NodeDataTypes>::type>,
					  						 public asc_NodeDown<KeyType, N, N, NodeDataTypes, EdgeDataTypes>
	{
		asc_Node(int id) : asc_NodeBase(id) {}
       
        friend std::ostream& operator<<(std::ostream& output, const asc_Node& node){
            output  << "Node(level=" << N
                    << ", id=" << node._node;
            if(node._down.size() > 0)
                for(auto it=node._down.cbegin(); it!=node._down.cend(); ++it)
                    output  << ", NodeDownID={'"
                            << it->first << "', "
                            << it->second->_node << "}";
            output << ")";
            return output;
        }
	};

	/*
	 * These are iterator adapters. The use of boost libraries is indicated.
	 */
	template <typename Iter, typename Data>
	struct node_iterator : public std::iterator<std::bidirectional_iterator_tag, Data> {
	public:
		using super = std::iterator<std::bidirectional_iterator_tag, Data>;
		node_iterator() {}
		node_iterator(Iter j) : i(j) {}
		node_iterator& operator++() { ++i; return *this; }
		node_iterator operator++(int) { auto tmp = *this; ++(*this); return tmp; }
		node_iterator& operator--() { --i; return *this; }
		node_iterator operator--(int) { auto tmp = *this; --(*this); return tmp; }
		bool operator==(node_iterator j) const { return i == j.i; }
		bool operator!=(node_iterator j) const { return !(*this == j); }
		typename super::reference operator*() { return i->second; }
        typename super::pointer operator->() { return i->second; }
	protected:
		Iter i;
	};

	template <typename Iter, typename Data>
	inline node_iterator<Iter,Data> make_node_iterator(Iter j) { return node_iterator<Iter,Data>(j); }

	template <typename Iter, typename Data>
	struct node_data_iterator : public std::iterator<std::bidirectional_iterator_tag, Data> {
	public:
		using super = std::iterator<std::bidirectional_iterator_tag, Data>;
		node_data_iterator() {}
		node_data_iterator(Iter j) : i(j) {}
		node_data_iterator& operator++() { ++i; return *this; }
		node_data_iterator operator++(int) { auto tmp = *this; ++(*this); return tmp; }
		node_data_iterator& operator--() { --i; return *this; }
		node_data_iterator operator--(int) { auto tmp = *this; --(*this); return tmp; }
		bool operator==(node_data_iterator j) const { return i == j.i; }
		bool operator!=(node_data_iterator j) const { return !(*this == j); }
		typename super::reference operator*() { return i->second->_data; }
		typename super::pointer operator->() { return i->second->_data; }
	protected:
		Iter i;
	};

	template <typename Iter, typename Data>
	inline 
	node_data_iterator<Iter,Data> make_node_data_iterator(Iter j) { return node_data_iterator<Iter,Data>(j); }
} // end namespace


template <typename K, typename... Ts>
struct simplicial_complex_traits_default
{
	template <std::size_t k> using all_void = int;
	using KeyType = K;
	using NodeTypes = util::type_holder<Ts...>;
	using EdgeTypes = typename util::int_type_map<std::size_t,
									                util::type_holder,
									                typename std::make_index_sequence<sizeof...(Ts)-1>,
									                all_void>::type;
};


template <typename traits>
class simplicial_complex {
public:
	using KeyType = typename traits::KeyType;
	using NodeDataTypes = typename traits::NodeTypes;
	using EdgeDataTypes = typename traits::EdgeTypes;
	using type_this = simplicial_complex<traits>;
	static const auto numLevels = NodeDataTypes::size;
	static const auto topLevel = numLevels-1;
	using LevelIndex = typename std::make_index_sequence<numLevels>;

	template <std::size_t k> using Node = detail::asc_Node<KeyType,k,topLevel,NodeDataTypes,EdgeDataTypes>;
	template <std::size_t k> using NodeData = typename util::type_get<k,NodeDataTypes>::type;
	template <std::size_t k> using EdgeData = typename util::type_get<k,EdgeDataTypes>::type;
	template <std::size_t k> using NodePtr = Node<k>*;


	simplicial_complex()
		: node_count(0)
	{
		_root = create_node<0>();
		for(auto& x : level_count)
		{
			x = 0;
		}
	}

	~simplicial_complex()
	{
		size_t count;
		remove_recurse<0,0>::apply(this, &_root, &_root + 1, count);
	}

	template <size_t n>
	void insert(const KeyType (&s)[n])
	{
		insert_for<0,n,false>::apply(this, _root, s);
	}

	template <size_t n>
	void insert(const KeyType (&s)[n], const NodeData<n>& data)
	{
		Node<n>* rval = insert_for<0,n,true>::apply(this, _root, s);
		rval->_data = data;
	}

	template <size_t n>
	std::array<KeyType,n> get_name(NodePtr<n> id)
	{
		std::array<KeyType,n> s;

		int i = 0;
		for(auto curr : id->_down)
		{
			s[i++] = curr.first;
		}
		
		return std::move(s);
	}

	std::array<KeyType,0> get_name(NodePtr<0> id)
	{
		std::array<KeyType,0> name;		
		return std::move(name);
	}

	// Node
	template <size_t n>
	NodePtr<n> get_id_up(const std::array<KeyType,n>& s)
	{
		return get_recurse<0,n>::apply(this, s.data(), _root);
	}

	template <size_t i, size_t j>
	auto get_id_up(NodePtr<i> nid, const std::array<KeyType,j>& s)
	{
		return get_recurse<i,j>::apply(this, s.data(), nid);
	}

	template <size_t i>
	auto get_id_up(NodePtr<i> nid, KeyType s)
	{
		return get_recurse<i,1>::apply(this, &s, nid);
	}

	template <size_t i, size_t j>
	auto get_id_down(NodePtr<i> nid, const std::array<KeyType,j>& s)
	{
		return get_down_recurse<i,j>::apply(this, s.data(), nid);
	}

	template <size_t i>
	auto get_id_down(NodePtr<i> nid, KeyType s)
	{
		return get_down_recurse<i,1>::apply(this, &s, nid);
	}

	/**
	 * @brief      Get the node data by node key
	 *
	 * @param[in]  s     An array of keys
	 *
	 * @tparam     n     The level of simplicial complex
	 *
	 * @return     The associate
	 */
	template <size_t n>
	NodeData<n>& get(const KeyType (&s)[n])
	{
		return get_recurse<0,n>::apply(this, s, _root)->_data;
	}

	template <size_t n>
	const NodeData<n>& get(const KeyType (&s)[n]) const
	{
		return get_recurse<0,n>::apply(this, s, _root)->_data;
	}

	template <size_t i>
	NodeData<i+1>& get(NodePtr<i> nid, KeyType s)
	{
		return get_recurse<i,1>::apply(this, &s, nid)->_data;
	}

	template <size_t i>
	NodeData<i>& get(NodePtr<i> nid)
	{
		return nid->_data;
	}

	NodeData<0>& get()
	{
		return _root->_data;
	}

	const NodeData<0>& get() const
	{
		return _root->_data;
	}

	template <size_t n, class Inserter>
	void get_cover(NodePtr<n> id, Inserter pos)
	{
		for(auto curr : id->_up)
		{
			pos++ = curr.first;
		}
	}

	template <size_t n>
	auto get_cover(NodePtr<n> id)
	{
		std::vector<KeyType> rval;
		get_cover(id, std::back_inserter(rval));
		return std::move(rval);
	}

	// Edge
	template <size_t n>
	EdgeData<n>& get_edge_up(NodePtr<n> nid, KeyType a)
	{
		return nid->_up[a]->_edge_data[a];
	}

	template <size_t n>
	EdgeData<n-1>& get_edge_down(NodePtr<n> nid, KeyType a)
	{
		return nid->_edge_data[a];
	}


	template <size_t n>
	bool exists(const KeyType (&s)[n]) const
	{
		
		return get_recurse<0,n>::apply(this, s, _root) != 0;
	}

	template <std::size_t k>
	auto size() const
	{
		return std::get<k>(levels).size();
	}

	template <std::size_t k>
	auto get_level_id()
	{
		auto begin = std::get<k>(levels).begin();
		auto end = std::get<k>(levels).end();
		auto data_begin = detail::make_node_iterator<decltype(begin),NodePtr<k>>(begin);
		auto data_end = detail::make_node_iterator<decltype(end),NodePtr<k>>(end);
		return util::make_range(data_begin, data_end);
	}

	template <std::size_t k>
	auto get_level_id() const
	{
		auto begin = std::get<k>(levels).cbegin();
		auto end = std::get<k>(levels).cend();
		auto data_begin = detail::make_node_iterator<decltype(begin), const NodePtr<k>>(begin);
		auto data_end = detail::make_node_iterator<decltype(end), const NodePtr<k>>(end);
		return util::make_range(data_begin, data_end);
	}

	template <std::size_t k>
	auto get_level()
	{
		auto begin = std::get<k>(levels).begin();
		auto end = std::get<k>(levels).end();
		auto data_begin = detail::make_node_data_iterator<decltype(begin),NodeData<k>>(begin);
		auto data_end = detail::make_node_data_iterator<decltype(end),NodeData<k>>(end);
		return util::make_range(data_begin, data_end);
	}

	template <std::size_t k>
	auto get_level() const
	{
		auto begin = std::get<k>(levels).cbegin();
		auto end = std::get<k>(levels).cend();
		auto data_begin = detail::make_node_data_iterator<decltype(begin), const NodeData<k>>(begin);
		auto data_end = detail::make_node_data_iterator<decltype(end), const NodeData<k>>(end);
		return util::make_range(data_begin, data_end);
	}

	template <size_t n>
	size_t remove(const KeyType (&s)[n])
	{
		Node<n>* root = get_recurse<0,n>::apply(this, s, _root);
		size_t count = 0;
		return remove_recurse<n,0>::apply(this, &root, &root + 1, count);
	}

	template <std::size_t k>
	void print_id()
	{
        std::cout << "level<" << k << ">.size()=" << this->size<k>() << std::endl; 
		
        auto ids = this->get_level_id<k>();
        for(auto& id : ids){
            std::cout << *id << std::endl;
        }
	}


private:
	/**
	 * Recursively deletes dependent nodes.
	 *
	 * @tparam     level  { description }
	 * @tparam     foo    { description }
	 */
	// The foo argument is necessary to get the compiler to shut up...
	template <size_t level, size_t foo>
	struct remove_recurse
	{
		template <typename T>
		static size_t apply(type_this* that, T begin, T end, size_t& count)
		{
			std::set<Node<level+1>*> next;
			// for each node of interest...
			for(auto i = begin; i != end; ++i)
			{
				auto up = (*i)->_up;
				for(auto j = up.begin(); j != up.end(); ++j)
				{
					next.insert(j->second);
				}
				that->remove_node(*i);
				++count;
			}
			return remove_recurse<level+1,foo>::apply(that, next.begin(), next.end(), count);
		}
	};

    // Terminal condition for remove_recurse
	template <size_t foo>
	struct remove_recurse<numLevels-1,foo>
	{
		template <typename T>
		static size_t apply(type_this* that, T begin, T end, size_t& count)
		{
			for(auto i = begin; i != end; ++i)
			{
				that->remove_node(*i);
				++count;
			}
			return count;
		}
	};

	template <size_t i, size_t n>
	struct get_recurse
	{
		static Node<i+n>* apply(const type_this* that, const KeyType* s, Node<i>* root)
		{
			if(root)
			{
				auto p = root->_up.find(*s);
				if(p != root->_up.end())
				{
					return get_recurse<i+1,n-1>::apply(that, s+1, root->_up[*s]);
				}
				else
				{
					return nullptr;
				}
			}
			else
			{
				return nullptr;
			}
		}
	};

	template <size_t i>
	struct  get_recurse<i,0>
	{
		static Node<i>* apply(const type_this* that, const KeyType* s, Node<i>* root)
		{
			return root;
		}
	};

	template <size_t i, size_t n>
	struct get_down_recurse
	{
		static Node<i-n>* apply(const type_this* that, const KeyType* s, Node<i>* root)
		{
			if(root)
			{
				auto p = root->_down.find(*s);
				if(p != root->_down.end())
				{
					return get_recurse<i-1,n-1>::apply(that, s+1, root->_down[*s]);
				}
				else
				{
					return nullptr;
				}
			}
			else
			{
				return nullptr;
			}
		}
	};

	template <size_t i>
	struct  get_down_recurse<i,0>
	{
		static Node<i>* apply(const type_this* that, const KeyType* s, Node<i>* root)
		{
			return root;
		}
	};


	template <size_t i, size_t n, bool do_insert>
	struct insert_for
	{
		static Node<i+n>* apply(type_this* that, Node<0>* root, const KeyType* begin)
		{
			that->extend<0,i,do_insert>(root, begin, *(begin+i));
			return insert_for<i+1,n-1,do_insert>::apply(that, root, begin);
		}
	};

	template <size_t i, bool do_insert>
	struct insert_for<i, 1, do_insert>
	{
		static Node<i+1>* apply(type_this* that, Node<0>* root, const KeyType* begin)
		{
			return that->extend<0,i,do_insert>(root, begin, *(begin+i));
		}
	};

	template <size_t level, size_t i, size_t n, bool do_insert>
	struct extend_for
	{
		static Node<level+i+n+1>* apply(type_this* that, Node<level>* root, const KeyType* begin, int value)
		{
			that->extend<level+1,i,do_insert>(root->_up[*(begin+i)], begin, value);
			return extend_for<level,i+1,n-1,do_insert>::apply(that, root, begin, value);
		}
	};

	template <size_t level, size_t i, bool do_insert>
	struct extend_for<level,i,1,do_insert>
	{
		static Node<level+i+1+1>* apply(type_this* that, Node<level>* root, const KeyType* begin, int value)
		{
			return that->extend<level+1,i,do_insert>(root->_up[*(begin+i)], begin, value);
		}
	};

	template <size_t level, size_t i, bool do_insert>
	struct extend_for<level,i,0,do_insert>
	{
		static Node<level+i+1>* apply(type_this* that, Node<level>* root, const KeyType* begin, int value)
		{
			return root->_up[value];
		}
	};

    /**
     *  @brief Backfill in the pointers from prior nodes to the new node
     *  @param root is a parent node
     *  @param nn is the new child node
     *  @param value is the exposed id of nn
     *  @return void
     */
	template <size_t level>
	void backfill(Node<level>* root, Node<level+1>* nn, int value)
	{
		for(auto curr = root->_down.begin(); curr != root->_down.end(); ++curr)
		{
			int v = curr->first;

			Node<level-1>* parent = curr->second;
			Node<level>* child = parent->_up[value]; 

			nn->_down[v] = child;
			child->_up[v] = nn;
		}
	}

    /**
     *  @brief Fill in the pointers from level 1 to 0.
     *  @param root is a level 0 node
     *  @param nn is a level 1 node
     *  @param value is the exposed id of nn
     *  @return void
     */
	void backfill(Node<0>* root, Node<1>* nn, int value)
	{
		return;
	}

	template <size_t level, size_t n, bool do_insert>
	Node<level+n+1>* extend(Node<level>* root, const KeyType* begin, int value)
	{
		if(root->_up.find(value) == root->_up.end())
		{
			Node<level+1>* nn;
			if(n == 0) // we are inserting the node the user requested.
			{
				nn = create_node<level+1>();
			}
			else // we are backfilling. The user may or may not want this.
			{
				nn = create_node<level+1>();
			}
			nn->_down[value] = root;
			root->_up[value] = nn;

			backfill(root, nn, value);
		}

		return extend_for<level,0,n,do_insert>::apply(this, root, begin, value);
	}

	template <size_t level>
	Node<level>* create_node()
	{
		auto p = new Node<level>(node_count++);
		++(level_count[level]);
	    
	    bool ret = std::get<level>(levels).insert(std::pair<size_t,NodePtr<level>>(node_count-1, p)).second; // node_count-1 to match the id's correctly
        // sanity check to make sure there aren't duplicate keys... 
        if (ret==false) {
            std::cout << "Error: Node '" << node_count << "' already existed with value " << *p << std::endl;
        }
		return p;
	}

	template <size_t level>
	void remove_node(Node<level>* p)
	{
		for(auto curr = p->_down.begin(); curr != p->_down.end(); ++curr)
		{
			curr->second->_up.erase(curr->first);
		}
		for(auto curr = p->_up.begin(); curr != p->_up.end(); ++curr)
		{
			curr->second->_down.erase(curr->first);
		}
		--(level_count[level]);
		std::get<level>(levels).erase(p->_node);
		delete p;
	}

	void remove_node(Node<0>* p)
	{
		for(auto curr = p->_up.begin(); curr != p->_up.end(); ++curr)
		{
			curr->second->_down.erase(curr->first);
		}
		--(level_count[0]);
		std::get<0>(levels).erase(p->_node);
		delete p;
	}

	void remove_node(Node<topLevel>* p)
	{
		for(auto curr = p->_down.begin(); curr != p->_down.end(); ++curr)
		{
			curr->second->_up.erase(curr->first);
		}
		--(level_count[topLevel]);
		std::get<topLevel>(levels).erase(p->_node);
		delete p;
	}

	Node<0>* _root;
	size_t node_count;
	std::array<size_t,numLevels> level_count;
	using NodePtrLevel = typename util::int_type_map<std::size_t, std::tuple, LevelIndex, NodePtr>::type;
	typename util::type_map<NodePtrLevel, detail::map>::type levels;
};

template <typename KeyType, typename... Ts>
using AbstractSimplicialComplex = simplicial_complex<simplicial_complex_traits_default<KeyType,Ts...>>;



template <class Complex, std::size_t level, class InsertIter>
void neighbors(Complex &F, typename Complex::template NodePtr<level> nid, InsertIter iter)
{
	for (auto a : F.get_name(nid))
	{
		auto id = F.get_id_down(nid,a);
		for(auto b : F.get_cover(id))
		{
			auto nbor = F.get_id_up(id,b);
			if(nbor != nid)
			{
				*iter++ = nbor;
			}
		}
	}
}

template <class Complex, std::size_t level, class InsertIter>
void neighbors_up(Complex &F, typename Complex::template NodePtr<level> nid, InsertIter iter)
{
	for (auto a : F.get_cover(nid))
	{
		auto id = F.get_id_up(nid,a);
		for(auto b : F.get_name(id))
		{
			auto nbor = F.get_id_down(id,b);
			if(nbor != nid)
			{
				*iter++ = nbor;
			}
		}
	}
}
