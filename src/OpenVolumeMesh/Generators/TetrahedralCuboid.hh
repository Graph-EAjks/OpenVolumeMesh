#pragma once

namespace OpenVolumeMesh {

struct SortedFace
{
    explicit SortedFace(std::vector<OpenVolumeMesh::VertexHandle> const& face)
        : v(face)
    {
        std::sort(v.begin(), v.end());
    }

    SortedFace(OpenVolumeMesh::VertexHandle v1,
            OpenVolumeMesh::VertexHandle v2,
            OpenVolumeMesh::VertexHandle v3)
        : v(3)
    {
        v[0] = v1;
        v[1] = v2;
        v[2] = v3;
        std::sort(v.begin(), v.end());
    }

    std::vector<OpenVolumeMesh::VertexHandle> v;
};

inline bool operator<(SortedFace const& f1, SortedFace const& f2)
{
    return std::lexicographical_compare(f1.v.begin(), f1.v.end(),
            f2.v.begin(), f2.v.end());
}

class TetrahedralCuboidGenerator
{
public:
    TetrahedralCuboidGenerator(PolyhedralMesh& mesh, Vector const& position, Vector const& length,
                                    unsigned const n_x, unsigned const n_y, unsigned const n_z);

private:
    void add_vertices(Vector const& position, Vector const& length);
    void get_cube_vertices(std::size_t i, std::size_t j, std::size_t k,
            std::vector<OpenVolumeMesh::VertexHandle>& v) const;

    void add_faces();
    void add_cube_type_1_faces(std::size_t i, std::size_t j, std::size_t k,
            std::vector<OpenVolumeMesh::VertexHandle> const& v);
    void add_cube_type_2_faces(std::size_t i, std::size_t j, std::size_t k,
            std::vector<OpenVolumeMesh::VertexHandle> const& v);

    void add_cells();
    void add_cube_type_1_cells(std::size_t i, std::size_t j, std::size_t k,
            std::vector<OpenVolumeMesh::VertexHandle> const& v);
    void add_cube_type_2_cells(std::size_t i, std::size_t j, std::size_t k,
            std::vector<OpenVolumeMesh::VertexHandle> const& v);

    PolyhedralMesh* mesh_;

    std::size_t size_[3];
    std::vector<OpenVolumeMesh::VertexHandle> vertices_;
    std::map<SortedFace, OpenVolumeMesh::FaceHandle> faces_;
};

template<typename MeshT>
FaceHandle add_axis_aligned_tetrahedral_cuboid(
        MeshT &_mesh,
        typename MeshT::Point const & _position,
        double _length,
        Vec3ui _counts)

} // namespace OpenVolumeMesh

