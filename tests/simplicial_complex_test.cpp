#include <iostream>
#include <memory>
#include <map>
#include <array>
#include <vector>
#include <queue>
#include <set>
#include <tuple>
#include <initializer_list>
#include <math.h>
#include <algorithm>
#include <tuple>
#include <future>
#include <experimental/optional>
#include "SimplicialComplex.h"


struct Vec3 {
	double x;
	double y;
	double z;
};

//AbstractSimplicialComplex* readOFF(const char *input_name);


//using ASC = AbstractSimplicialComplex<size_t,int,Vec3,int,int>;
//ASC* readOFF(const char *input_name);
struct Orientable {
	int orientation;
};

struct Face : Orientable {
	Face(int orient, long m) : Orientable{orient}, mark(m) {}
	long mark;
};

struct complex_traits
{
	using KeyType = size_t;
	using NodeTypes = util::type_holder<int,Vec3,int,Orientable>;
	using EdgeTypes = util::type_holder<Orientable,Orientable,Orientable>;
};
using ASC = simplicial_complex<complex_traits>;
ASC* readOFF(const char *input_name);

void print_vertex(const ASC& F)
{
	for(auto curr : F.get_level<1>())
	{
		std::cout << curr.x << " " << curr.y << " " << curr.z << std::endl;
	}
}

template <std::size_t level, class InsertIter>
void neighbors(ASC &F, ASC::NodeId<level> nid, InsertIter iter)
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

template <std::size_t level, class InsertIter>
void neighbors_up(ASC &F, ASC::NodeId<level> nid, InsertIter iter)
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

template <class Complex, class SizeT >
struct init_orientation_helper {};

template <class Complex, std::size_t k >
struct init_orientation_helper<Complex, std::integral_constant<std::size_t, k>>
{
	static void f(Complex& F)
	{
		for(auto curr : F.template get_level_id<k>())
		{
			for(auto a : F.get_cover(curr))
			{
				int orient = 1;
				for(auto b : F.get_name(curr))
				{
					if(a > b)
					{
						orient *= -1;
					}
					else
					{
						break;
					}
				}
				F.get_edge_up(curr,a).orientation = orient;
			}
		}

		init_orientation_helper<Complex,std::integral_constant<std::size_t,k+1>>::f(F);
	}
};

template <typename Complex>
struct init_orientation_helper<Complex, std::integral_constant<std::size_t, Complex::topLevel>>
{
	static void f(Complex& F) {}
};

template <typename Complex>
void init_orientation(Complex& F)
{
	init_orientation_helper<Complex,std::integral_constant<std::size_t,0>>::f(F);
}

template <typename Complex>
void clear_orientation(Complex& F)
{
	constexpr std::size_t k = Complex::topLevel;

	for(auto& curr : F.template get_level<k>())
	{
		curr.orientation = 0;
	}
}

template <typename Complex>
void compute_orientation(Complex& F)
{
	constexpr std::size_t k = Complex::topLevel - 1;

	std::queue<ASC::NodeId<2>> frontier;
	std::set<ASC::NodeId<2>> visited;
	int connected_components = 0;
	bool orientable = true;
	bool psuedo_manifold = true;
	for(auto outer : F.template get_level_id<k>())
	{
//		if(!F.get(outer))
		if(visited.find(outer) == visited.end())
		{
			++connected_components;
			frontier.push(outer);

			while(!frontier.empty())
			{
				auto curr = frontier.front();
//				if(!F.get(curr))
				if(visited.find(curr) == visited.end())
				{
//					F.get(curr) = 1;
					visited.insert(curr);

					auto w = F.get_cover(curr);

					if(w.size() == 1)
					{
						std::cout << curr << ":" << w[0] << " ~ Boundary" << std::endl;
					}
					else if(w.size() == 2)
					{
						auto& edge0 = F.get_edge_up(curr, w[0]);
						auto& edge1 = F.get_edge_up(curr, w[1]);

						auto& node0 = F.get(curr, w[0]);
						auto& node1 = F.get(curr, w[1]);

						if(node0.orientation == 0)
						{
							if(node1.orientation == 0)
							{
								node0.orientation = 1;
								node1.orientation = -edge1.orientation * edge0.orientation * node0.orientation;
							}
							else
							{
								node0.orientation = -edge0.orientation * edge1.orientation * node1.orientation;
							}
						}
						else
						{
							if(node1.orientation == 0)
							{
								node1.orientation = -edge1.orientation * edge0.orientation * node0.orientation;
							}
							else
							{
								if(edge0.orientation*node0.orientation + edge1.orientation*node1.orientation != 0)
								{
									orientable = false;
									std::cout << "+++++" << std::endl;
									std::cout << edge0.orientation << " : " << node0.orientation << std::endl;
									std::cout << edge1.orientation << " : " << node1.orientation << std::endl;

									std::cout << " : "
									          << edge0.orientation*node0.orientation + edge1.orientation*node1.orientation
									          << std::endl;
									std::cout << "-----"
									          << std::endl;
									std::cout << "Non-Orientable: "
									          << edge0.orientation*node0.orientation + edge1.orientation*node1.orientation
									          << std::endl;
								}
							}
						}

						std::vector<ASC::NodeId<2>> tmp;
						neighbors_up(F, curr, std::back_inserter(tmp));
						for(auto nid : tmp)
						{
							frontier.push(nid);
						}
					}
					else
					{
						psuedo_manifold = false;
					}
				}
				frontier.pop();
			}
		}
	}

	std::cout << "Connected Components: " << connected_components << std::endl;
	std::cout << "Orientable: " << orientable << std::endl;
	std::cout << "Psuedo-manifold: " << psuedo_manifold << std::endl;
}




int main(int argc, char* argv[])
{
	if(argc != 2)
	{
		std::cerr << "Wrong arguments passed" << std::endl;
		return -1;
	}
	auto pF = readOFF(argv[1]);

	init_orientation(*pF);

	clear_orientation(*pF);
	compute_orientation(*pF);


	std::map<ASC::NodeId<3>,int> X;
	std::map<ASC::NodeId<2>,int> Y;
	std::map<ASC::NodeId<1>,int> Z;

	for(auto curr : pF->get_level_id<3>())
	{
		X[curr] = pF->get(curr).orientation;
	}

	for(auto curr : pF->get_level_id<3>())
	{
		std::array<ASC::KeyType,3> w = pF->get_name(curr);
		for(auto a : w)
		{
			int A = pF->get_edge_down(curr, a).orientation;
			ASC::NodeId<2> id = pF->get_id_down(curr, a);

			if(Y.find(id) != Y.end())
			{
				Y[id] += A*X[curr];
			}
			else
			{
				Y[id] = A*X[curr];
			}
		}
	}

	for(auto curr : pF->get_level_id<2>())
	{
		std::array<ASC::KeyType,2> w = pF->get_name(curr);
		for(auto a : w)
		{
			int A = pF->get_edge_down(curr, a).orientation;
			ASC::NodeId<1> id = pF->get_id_down(curr, a);

			if(Z.find(id) != Z.end())
			{
				Z[id] += A*Y[curr];
			}
			else
			{
				Z[id] = A*Y[curr];
			}
		}
	}

	for(auto curr : pF->get_level_id<2>())
	{
		if(Y[curr] != 0)
		{
			std::cout << curr << " : " << Y[curr] << std::endl;
		}
	}
	/*
	count = 0;
	for(auto curr : pF->get_level_id<1>())
	{
		std::cout << Z[curr] << std::endl;
	}
	*/
	delete pF;

	return -1;
}

int rgb_to_marker(float r, float g, float b)
{
    if (r < 0 || r > 1 || g < 0 || g > 1 || b < 0 || b > 1)
    {
        printf("Expected individual RGB value to be betwen 0 and 1.\n");
        exit(1);
    }

    return round(r * 10) * 121 + round(g * 10) * 11 + round(b * 10);
}


ASC* readOFF(const char *input_name)
{
	ASC* F = new ASC();
    unsigned int n, m;
    unsigned int a, b, c;
    float x, y, z, color_r, color_g, color_b, color_a;
    unsigned int          v_n, t_n, e_n, character;
    char         line[256];
    FILE        *fin = NULL;
    fpos_t       fp;
    bool         has_marker = false;
    std::map<Simplex, Vec3> vertex_position;
    std::map<Simplex, int> face_marking;


    if ((fin = fopen(input_name, "r")) == NULL)
    {
        std::cerr << "Read error. File \'" << input_name << "\' could not be read." << std::endl;
        return NULL;
    }

    while (fgets(line, 256, fin) != NULL)
    {
        if ((line[0] == 'O') && (line[1] == 'F') && (line[2] == 'F'))
        {
            break;
        }
    }

    if (fscanf(fin, "%d %d %d\n", &v_n, &t_n, &e_n) != 3)
    {

        std::cerr << "Read error. Expected 3 integer number for the number of vertices and simplices." << std::endl;
        return NULL;
    }

    std::cout << "   vertices: " << v_n << " --- simplices: " << t_n << std::endl;

    // Read vertex coordinates
    for (n = 0; n < v_n; n++) // surfmesh->num_vertices; n++) {
    {
        if (fscanf(fin, "%f %f %f\n", &x, &y, &z) != 3)
        {
            std::cerr << "Read error. Expected 3 floats for the coordinates of vertex:" << n << std::endl;
            return NULL;
        }

        //        surfmesh->vertex.push_back(FLTVECT(x,y,z));

        Vec3 v{x,y,z};
        F->insert<1>({n}, v);
    }

    // Read the number of vertices per simplex from the first line of simplices
    n = fscanf(fin, "%d", &m);

    std::cout << m << std::endl;

    // Check format of the read simplices
    if ((m != 3) && (m != 4))
    {
        std::cerr << "Read error. Expected a 3 or 4 for the first value in the first simplex line." << std::endl;
        return NULL;
    }

    // Input is surface mesh
    else if (m == 3)
    {
        std::cout << "   Input is surface mesh." << std::endl;

        // Read the rest of the first line
        if (fscanf(fin, "%d %d %d", &a, &b, &c) != 3)
        {
            std::cerr << "Read error. Expected 3 integers for the first simplex." << std::endl;
            return NULL;
        }

        // Get position of file and check if we have more symbols
        fgetpos(fin, &fp);

        // Read any blanks
        character = fgetc(fin);

        while (character == ' ')
        {
            character = fgetc(fin);
        }

        // If we do not have en end of line character we have marker
        if (character != '\n')
        {
            // Set marker flag
            has_marker = true;

            // Rewind the file pointer
            fsetpos(fin, &fp);

            // Read first rgb values
            if (fscanf(fin, "%f %f %f %f", &color_r, &color_g, &color_b, &color_a) != 4)
            {
                std::cerr << "Read error. Expected 4 floats for the RGBA values of the first face." << std::endl;
                return NULL;
            }

            while (fgetc(fin) != '\n')
            {}

        	F->insert<3>({a,b,c}, Face(0, rgb_to_marker(color_r, color_g, color_b)));
        }
        else
        {
	        F->insert<3>({a,b,c});
        }

        // Read the rest of the simplices
        for (n = 1; n < t_n; n++)
        {
            if (fscanf(fin, "%d %d %d %d", &m, &a, &b, &c) != 4)
            {
                std::cerr << "Read error. Expected 4 integers for simplex " << n << std::endl;
                return nullptr;
            }

            // If we have face markers
            if (has_marker)
            {
                // Read first rgb values
                if (fscanf(fin, "%f %f %f %f", &color_r, &color_g, &color_b,
                           &color_a) != 4)
                {
                    std::cerr << "Read error. Expected 4 floats for the RGBA values of face:" << n << std::endl;
                    return NULL;
                }

                // Get the other face markers
                //surfmesh->face[n].m = rgb_to_marker(color_r, color_g, color_b);
	        	F->insert<3>({a,b,c}, Face(0, rgb_to_marker(color_r, color_g, color_b)));
            }
            else
            {
		        F->insert<3>({a,b,c});
            }

            // Skip any additional character on this line
            while (fgetc(fin) != '\n')
            {}
        }

    }

    // Close file
    fclose(fin);
    /*
    // Input is tetrahedral mesh.
    else
    {
        printf("   Input is volumetric mesh...\n");

        // Read the rest of the first line
        if (fscanf(fin, "%d %d %d %d", &a, &b, &c, &d) != 4)
        {
            printf("Read error. Expected 4 integers for the first simplex.\n");
            return nullptr;
        }

        // Skip any additional character on this line
        while (fgetc(fin) != '\n')
        {}

        // Assign first simplex from allready read data
    	F->insert({a,b,c,d});

        // Read the rest of the simplices
        for (n = 1; n < t_n; n++)
        {
            if (fscanf(fin, "%d %d %d %d %d", &m, &a, &b, &c, &d) != 5)
            {
                printf("Read error. Expected 5 integers for simplex %d.\n", n);
                return nullptr;
            }
	    	F->insert({a,b,c,d});

            // Skip any additional character on this line
            while (fgetc(fin) != '\n')
            {}
        }

        fclose(fin);
    }
    */

    return F;
}
