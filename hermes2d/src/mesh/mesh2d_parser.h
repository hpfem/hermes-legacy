// mesh2d_parser.h

# include <iostream>
# include <string>
# include <fstream>
# include <sstream>
# include <vector>
# include <cstdlib>
# include <cassert>
# include <map>

class MeshData
{
	std::string mesh_file_;
	void strip(std::string &str);
	
public:
	std::map< std::string, std::vector< std::string > > vars_;

	int n_vert, n_el, n_bdy, n_curv, n_ref;

	std::vector<double> x_vertex, y_vertex;
	
	std::vector<int> en1, en2, en3, en4;
	std::vector<std::string> e_mtl;
	
	std::vector<int> bdy_first, bdy_second;
	std::vector<std::string> bdy_type;
	
	std::vector<int> curv_first, curv_second;
	std::vector<double> curv_third;
	std::vector<std::string> curv_inner_pts;
	std::vector<std::string> curv_knots;
	std::vector<bool> curv_nurbs;
	
	std::vector<int> ref_elt;
	std::vector<int> ref_type;
	
	void parse_mesh(void);
	
	MeshData(const std::string &mesh_file) 
	{
		mesh_file_ = mesh_file;
	}
	
	MeshData(const MeshData &m) : mesh_file_(m.mesh_file_), n_vert(m.n_vert), n_el(m.n_el), n_bdy(m.n_bdy), n_curv(m.n_curv), n_ref(m.n_ref)
	{
		vars_ = m.vars_;

		x_vertex = m.x_vertex;
		y_vertex = m.y_vertex;
		
		en1 = m.en1; en2 = m.en2; en3 = m.en3; en4 = m.en4;
		e_mtl = m.e_mtl;
		
		bdy_first = m.bdy_first; bdy_second = m.bdy_second;
		bdy_type = m.bdy_type;
		
		curv_first = m.curv_first; curv_second = m.curv_second;
		curv_third = m.curv_third;
		curv_inner_pts = m.curv_inner_pts;
		curv_knots = m.curv_knots;
		curv_nurbs = m.curv_nurbs;
		
		ref_elt = m.ref_elt;
		ref_type = m.ref_type;
	}
	
	MeshData& operator = (const MeshData &m)
	{
		assert(&m != this);
		
		mesh_file_ = m.mesh_file_;
		vars_ = m.vars_;
		
		n_vert = m.n_vert;
		n_el = m.n_el;
		n_bdy = m.n_bdy;
		n_curv = m.n_curv;
		n_ref = m.n_ref;
		
		x_vertex = m.x_vertex;
		y_vertex = m.y_vertex;
		
		en1 = m.en1; en2 = m.en2; en3 = m.en3; en4 = m.en4;
		e_mtl = m.e_mtl;
		
		bdy_first = m.bdy_first; bdy_second = m.bdy_second;
		bdy_type = m.bdy_type;
		
		curv_first = m.curv_first; curv_second = m.curv_second;
		curv_third = m.curv_third;
		curv_inner_pts = m.curv_inner_pts;
		curv_knots = m.curv_knots;
		curv_nurbs = m.curv_nurbs;
		
		ref_elt = m.ref_elt;
		ref_type = m.ref_type;
		
		return *this;
	}
	
	~MeshData() {}
};
