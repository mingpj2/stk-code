#ifndef STKINSTANCEDSCENENODE_HPP
#define STKINSTANCEDSCENENODE_HPP

#include "stkmesh.hpp"
#include "utils/leak_check.hpp"

class ListInstancedMatDefault : public MeshList<ListInstancedMatDefault, GLMesh *, size_t>
{};

class ListInstancedMatAlphaRef : public MeshList<ListInstancedMatAlphaRef, GLMesh *, size_t>
{};

class ListInstancedMatGrass : public MeshList<ListInstancedMatGrass, GLMesh *, size_t, core::vector3df, core::vector3df>
{};

class ListInstancedMatNormalMap : public MeshList<ListInstancedMatNormalMap, GLMesh *, size_t>
{};

class STKInstancedSceneNode : public irr::scene::CMeshSceneNode
{
protected:
    int m_ref_count;
    std::vector<GLMesh *> MeshSolidMaterial[MAT_COUNT];
    std::vector<GLMesh> GLmeshes;
    std::vector<float> instance_pos;
    core::matrix4 ModelViewProjectionMatrix, TransposeInverseModelView;
    GLuint instances_vbo;
    void createGLMeshes();
    bool isMaterialInitialized;
    void setFirstTimeMaterial();
    void initinstancedvaostate(GLMesh &mesh);
    void cleanGL();
    core::vector3df windDir;
public:
    STKInstancedSceneNode(irr::scene::IMesh* mesh, ISceneNode* parent, irr::scene::ISceneManager* mgr, irr::s32 id,
        const irr::core::vector3df& position = irr::core::vector3df(0, 0, 0),
        const irr::core::vector3df& rotation = irr::core::vector3df(0, 0, 0),
        const irr::core::vector3df& scale = irr::core::vector3df(1.0f, 1.0f, 1.0f));
    ~STKInstancedSceneNode();
    virtual void render();
    void addInstance(const core::vector3df &origin, const core::vector3df &orientation, const core::vector3df &scale);

    int getInstanceCount() const { return instance_pos.size() / 9; }

    core::matrix4 getInstanceTransform(int id);

    void instanceGrab() { m_ref_count++; }
    void instanceDrop()
    {
        m_ref_count--;
        if (m_ref_count <= 0)
        {
            delete this;
        }
    }

    LEAK_CHECK();
};

#endif
