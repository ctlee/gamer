#pragma once

#include <iostream>
#include <fstream>
#include <ostream>
#include "SimplicialComplex.h"
#include "SimplicialComplexVisitors.h"

/**
 * @brief      Returns a string representation of the vertex subsimplicies of a given simplex
 *
 * @param[in]  A     Array containing name of a simplex
 *
 * @tparam     T     Type of the array elements
 * @tparam     k     Number of elements
 *
 * @return     String representation of the object.
 */
template <typename T, std::size_t k>
std::string to_string(const std::array<T,k>& A)
{
    if (k==0){
        return "{root}";
    }
	std::string out;
    out += "\"{";
    for(int i = 0; i + 1 < k; ++i)
    {
        out += std::to_string(A[i]) + ",";
    }
    if(k > 0)
    {
        out += std::to_string(A[k-1]);
    }
    out += "}\"";
    return out;
}

/**
 * @brief      Visitor for printing connectivity of the ASC
 *
 * @tparam     Complex  Type details of the ASC
 */
template <typename Complex>
struct GraphVisitor
{
    std::ostream& fout;
    GraphVisitor(std::ostream& os) : fout(os) {}

    template <std::size_t level> 
    bool visit(const Complex& F, typename Complex::template SimplexID<level> s)
    {
    	auto name = to_string(F.get_name(s));

   		auto covers = F.get_cover(s);
   		for(auto cover : covers){
   			auto edge = F.get_edge_up(s, cover);
   			auto nextName = to_string(F.get_name(edge.up()));
   			if((*edge).orientation == 1){
   				fout    << "   " << name << " -> " 
   					    << nextName << std::endl;
   			}
   			else{
				fout    << "   " << nextName << " -> " 
				   		<< name << std::endl;
   			}
	    }
    	return true;
    }

    bool visit(const Complex& F, typename Complex::template SimplexID<Complex::topLevel-1> s){

        auto name = to_string(F.get_name(s));
        auto covers = F.get_cover(s);
        for(auto cover : covers){
            auto edge = F.get_edge_up(s, cover);
            auto nextName = to_string(F.get_name(edge.up()));
            auto orient = (*edge.up()).orientation;
            if (orient == 1){
                nextName = "+ " + nextName;
            }
            else{
                nextName = "- " + nextName;
            }
            if((*edge).orientation == 1){
                fout    << "   " << name << " -> " 
                        << nextName << std::endl;
            }
            else{
                fout    << "   " << nextName << " -> " 
                        << name << std::endl;
            }
        }
        return true;
    }

    void visit(const Complex& F, typename Complex::template SimplexID<Complex::topLevel> s){}
};

template <typename Complex, typename K>
struct DotHelper {};

/**
 * @brief      Print out the names of simplices at each level
 *
 * @tparam     Complex  Type of the ASC
 * @tparam     k        Level to traverse
 */
template <typename Complex, std::size_t k>
struct DotHelper<Complex, std::integral_constant<std::size_t, k>>
{
    static void printlevel(std::ofstream& fout, const Complex& F){
        auto nodes = F.template get_level_id<k>();
        fout    << "subgraph cluster_" << k << " {\n"
                << "label=\"Level " << k << "\"\n";
        for (auto node : nodes){
            fout << to_string(F.get_name(node)) << ";";
        }
        fout    << "\n}\n";
        DotHelper<Complex, std::integral_constant<std::size_t,k+1>>::printlevel(fout, F);
    }
};

/**
 * @brief      Print out the names of simplices at the top level
 *
 * @tparam     Complex  Type of the ASC
 */
template <typename Complex>
struct DotHelper<Complex, std::integral_constant<std::size_t, Complex::topLevel>>
{
    static void printlevel(std::ofstream& fout, const Complex& F){
        auto nodes = F.template get_level_id<Complex::topLevel>();
        fout    << "subgraph cluster_" << Complex::topLevel << " {\n"
                << "label=\"Level " << Complex::topLevel << "\"\n";
        for (auto node : nodes){
            fout << to_string(F.get_name(node)) << ";";
        }
        fout    << "\n}\n";
    }
};


/**
 * @brief      Writes out the topology of an ASC into the dot format. The resulting dot file can be rendered into an image using tools such as GraphViz. dot -Tpng file.dot > image.png.
 *
 * @param[in]  filename  The filename
 * @param[in]  F         ASC of interest
 *
 * @tparam     Complex   Type of the ASC
 */
template <typename Complex>
void writeDOT(const std::string& filename, const Complex& F){
	std::ofstream fout(filename);
	if(!fout.is_open())
	{
	    std::cerr   << "File '" << filename 
	                << "' could not be writen to." << std::endl;
        fout.close();
	    exit(1); 
	}

	fout 	<< "digraph {\n"
            << "node [shape = record,height = .1]"
            << "splines=line;\n"
			<< "dpi=300;\n";
    auto v = GraphVisitor<Complex>(fout);
	visit_node_up(v, F, F.get_simplex_up());

    DotHelper<Complex, std::integral_constant<std::size_t, 0>>::printlevel(fout, F);

	fout << "}\n";
    fout.close();
}