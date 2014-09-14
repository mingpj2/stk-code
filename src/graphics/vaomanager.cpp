#include "vaomanager.hpp"
#include "stkmesh.hpp"

VAOManager::VAOManager()
{
    vao[0] = vao[1] = vao[2] = 0;
    vbo[0] = vbo[1] = vbo[2] = 0;
    ibo[0] = ibo[1] = ibo[2] = 0;
    last_vertex[0] = last_vertex[1] = last_vertex[2] = 0;
    last_index[0] = last_index[1] = last_index[2] = 0;
    RealVBOSize[0] = RealVBOSize[1] = RealVBOSize[2] = 0;
    RealIBOSize[0] = RealIBOSize[1] = RealIBOSize[2] = 0;

    for (unsigned i = 0; i < InstanceTypeCount; i++)
    {
        glGenBuffers(1, &instance_vbo[i]);
        glBindBuffer(GL_ARRAY_BUFFER, instance_vbo[i]);
        if (irr_driver->hasBufferStorageExtension())
        {
            glBufferStorage(GL_ARRAY_BUFFER, 10000 * sizeof(InstanceData), 0, GL_MAP_WRITE_BIT | GL_MAP_PERSISTENT_BIT);
            Ptr[i] = glMapBufferRange(GL_ARRAY_BUFFER, 0, 10000 * sizeof(InstanceData), GL_MAP_WRITE_BIT | GL_MAP_PERSISTENT_BIT);
        }
        else
        {
            glBufferData(GL_ARRAY_BUFFER, 10000 * sizeof(InstanceData), 0, GL_STREAM_DRAW);
        }
    }
}

void VAOManager::cleanInstanceVAOs()
{
    std::map<std::pair<video::E_VERTEX_TYPE, InstanceType>, GLuint>::iterator It = InstanceVAO.begin(), E = InstanceVAO.end();
    for (; It != E; It++)
        glDeleteVertexArrays(1, &(It->second));
    InstanceVAO.clear();
}

VAOManager::~VAOManager()
{
    cleanInstanceVAOs();
    for (unsigned i = 0; i < 3; i++)
    {
        if (vbo[i])
            glDeleteBuffers(1, &vbo[i]);
        if (ibo[i])
            glDeleteBuffers(1, &ibo[i]);
        if (vao[i])
            glDeleteVertexArrays(1, &vao[i]);
    }
    for (unsigned i = 0; i < InstanceTypeCount; i++)
    {
        glDeleteBuffers(1, &instance_vbo[i]);
    }

}

void VAOManager::regenerateBuffer(enum VTXTYPE tp, size_t newlastvertex, size_t newlastindex)
{
    glBindVertexArray(0);

    if (newlastindex >= RealVBOSize[tp])
    {
        while (newlastindex >= RealVBOSize[tp])
            RealVBOSize[tp] = 2 * RealVBOSize[tp] + 1;
        GLuint newVBO;
        glGenBuffers(1, &newVBO);
        glBindBuffer(GL_ARRAY_BUFFER, newVBO);
        if (irr_driver->hasBufferStorageExtension())
        {
            glBufferStorage(GL_ARRAY_BUFFER, RealVBOSize[tp] * getVertexPitch(tp), 0, GL_MAP_WRITE_BIT | GL_MAP_PERSISTENT_BIT);
            VBOPtr[tp] = glMapBufferRange(GL_ARRAY_BUFFER, 0, RealVBOSize[tp] * getVertexPitch(tp), GL_MAP_WRITE_BIT | GL_MAP_PERSISTENT_BIT);
        }
        else
            glBufferData(GL_ARRAY_BUFFER, RealVBOSize[tp] * getVertexPitch(tp), 0, GL_DYNAMIC_DRAW);

        if (vbo[tp])
        {
            GLuint oldVBO = vbo[tp];
            glBindBuffer(GL_COPY_WRITE_BUFFER, newVBO);
            glBindBuffer(GL_COPY_READ_BUFFER, oldVBO);
            glCopyBufferSubData(GL_COPY_READ_BUFFER, GL_COPY_WRITE_BUFFER, 0, 0, last_vertex[tp] * getVertexPitch(tp));
            glDeleteBuffers(1, &oldVBO);
        }
        vbo[tp] = newVBO;
    }
    last_vertex[tp] = newlastvertex;

    if (newlastindex >= RealIBOSize[tp])
    {
        while (newlastindex >= RealIBOSize[tp])
            RealIBOSize[tp] = 2 * RealIBOSize[tp] + 1;
        GLuint newIBO;
        glGenBuffers(1, &newIBO);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, newIBO);
        if (irr_driver->hasBufferStorageExtension())
        {
            glBufferStorage(GL_ELEMENT_ARRAY_BUFFER, RealIBOSize[tp] * sizeof(u16), 0, GL_MAP_WRITE_BIT | GL_MAP_PERSISTENT_BIT);
            IBOPtr[tp] = glMapBufferRange(GL_ELEMENT_ARRAY_BUFFER, 0, RealIBOSize[tp] * sizeof(u16), GL_MAP_WRITE_BIT | GL_MAP_PERSISTENT_BIT);
        }
        else
            glBufferData(GL_ELEMENT_ARRAY_BUFFER, RealIBOSize[tp] * sizeof(u16), 0, GL_STATIC_DRAW);

        if (ibo[tp])
        {
            GLuint oldIBO = ibo[tp];
            glBindBuffer(GL_COPY_WRITE_BUFFER, newIBO);
            glBindBuffer(GL_COPY_READ_BUFFER, oldIBO);
            glCopyBufferSubData(GL_COPY_READ_BUFFER, GL_COPY_WRITE_BUFFER, 0, 0, last_index[tp] * sizeof(u16));
            glDeleteBuffers(1, &oldIBO);
        }
        ibo[tp] = newIBO;
    }
    last_index[tp] = newlastindex;
}

void VAOManager::regenerateVAO(enum VTXTYPE tp)
{
    if (vao[tp])
        glDeleteVertexArrays(1, &vao[tp]);
    glGenVertexArrays(1, &vao[tp]);
    glBindVertexArray(vao[tp]);
    glBindBuffer(GL_ARRAY_BUFFER, vbo[tp]);
    switch (tp)
    {
    case VTXTYPE_STANDARD:
        // Position
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, getVertexPitch(tp), 0);
        // Normal
        glEnableVertexAttribArray(1);
        glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, getVertexPitch(tp), (GLvoid*)12);
        // Color
        glEnableVertexAttribArray(2);
        glVertexAttribPointer(2, 4, GL_UNSIGNED_BYTE, GL_TRUE, getVertexPitch(tp), (GLvoid*)24);
        // Texcoord
        glEnableVertexAttribArray(3);
        glVertexAttribPointer(3, 2, GL_FLOAT, GL_FALSE, getVertexPitch(tp), (GLvoid*)28);
        break;
    case VTXTYPE_TCOORD:
        // Position
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, getVertexPitch(tp), 0);
        // Normal
        glEnableVertexAttribArray(1);
        glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, getVertexPitch(tp), (GLvoid*)12);
        // Color
        glEnableVertexAttribArray(2);
        glVertexAttribPointer(2, 4, GL_UNSIGNED_BYTE, GL_TRUE, getVertexPitch(tp), (GLvoid*)24);
        // Texcoord
        glEnableVertexAttribArray(3);
        glVertexAttribPointer(3, 2, GL_FLOAT, GL_FALSE, getVertexPitch(tp), (GLvoid*)28);
        // SecondTexcoord
        glEnableVertexAttribArray(4);
        glVertexAttribPointer(4, 2, GL_FLOAT, GL_FALSE, getVertexPitch(tp), (GLvoid*)36);
        break;
    case VTXTYPE_TANGENT:
        // Position
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, getVertexPitch(tp), 0);
        // Normal
        glEnableVertexAttribArray(1);
        glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, getVertexPitch(tp), (GLvoid*)12);
        // Color
        glEnableVertexAttribArray(2);
        glVertexAttribPointer(2, 4, GL_UNSIGNED_BYTE, GL_TRUE, getVertexPitch(tp), (GLvoid*)24);
        // Texcoord
        glEnableVertexAttribArray(3);
        glVertexAttribPointer(3, 2, GL_FLOAT, GL_FALSE, getVertexPitch(tp), (GLvoid*)28);
        // Tangent
        glEnableVertexAttribArray(5);
        glVertexAttribPointer(5, 3, GL_FLOAT, GL_FALSE, getVertexPitch(tp), (GLvoid*)36);
        // Bitangent
        glEnableVertexAttribArray(6);
        glVertexAttribPointer(6, 3, GL_FLOAT, GL_FALSE, getVertexPitch(tp), (GLvoid*)48);
        break;
    }

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo[tp]);
    glBindVertexArray(0);
}

void VAOManager::regenerateInstancedVAO()
{
    cleanInstanceVAOs();

    enum video::E_VERTEX_TYPE IrrVT[] = { video::EVT_STANDARD, video::EVT_2TCOORDS, video::EVT_TANGENTS };
    for (unsigned i = 0; i < VTXTYPE_COUNT; i++)
    {
        video::E_VERTEX_TYPE tp = IrrVT[i];
        if (!vbo[tp] || !ibo[tp])
            continue;
        GLuint vao = createVAO(vbo[tp], ibo[tp], tp);
        glBindBuffer(GL_ARRAY_BUFFER, instance_vbo[InstanceTypeDefault]);

        glEnableVertexAttribArray(7);
        glVertexAttribPointer(7, 3, GL_FLOAT, GL_FALSE, sizeof(InstanceData), 0);
        glVertexAttribDivisor(7, 1);
        glEnableVertexAttribArray(8);
        glVertexAttribPointer(8, 3, GL_FLOAT, GL_FALSE, sizeof(InstanceData), (GLvoid*)(3 * sizeof(float)));
        glVertexAttribDivisor(8, 1);
        glEnableVertexAttribArray(9);
        glVertexAttribPointer(9, 3, GL_FLOAT, GL_FALSE, sizeof(InstanceData), (GLvoid*)(6 * sizeof(float)));
        glVertexAttribDivisor(9, 1);
        glEnableVertexAttribArray(10);
        glVertexAttribIPointer(10, 2, GL_UNSIGNED_INT, sizeof(InstanceData), (GLvoid*)(9 * sizeof(float)));
        glVertexAttribDivisor(10, 1);
        glEnableVertexAttribArray(11);
        glVertexAttribIPointer(11, 2, GL_UNSIGNED_INT, sizeof(InstanceData), (GLvoid*)(9 * sizeof(float) + 2 * sizeof(unsigned)));
        glVertexAttribDivisor(11, 1);
        InstanceVAO[std::pair<video::E_VERTEX_TYPE, InstanceType>(tp, InstanceTypeDefault)] = vao;

        vao = createVAO(vbo[tp], ibo[tp], tp);
        glBindBuffer(GL_ARRAY_BUFFER, instance_vbo[InstanceTypeShadow]);

        glEnableVertexAttribArray(7);
        glVertexAttribPointer(7, 3, GL_FLOAT, GL_FALSE, sizeof(InstanceData), 0);
        glVertexAttribDivisor(7, 1);
        glEnableVertexAttribArray(8);
        glVertexAttribPointer(8, 3, GL_FLOAT, GL_FALSE, sizeof(InstanceData), (GLvoid*)(3 * sizeof(float)));
        glVertexAttribDivisor(8, 1);
        glEnableVertexAttribArray(9);
        glVertexAttribPointer(9, 3, GL_FLOAT, GL_FALSE, sizeof(InstanceData), (GLvoid*)(6 * sizeof(float)));
        glVertexAttribDivisor(9, 1);
        glEnableVertexAttribArray(10);
        glVertexAttribIPointer(10, 2, GL_UNSIGNED_INT, sizeof(InstanceData), (GLvoid*)(9 * sizeof(float)));
        glVertexAttribDivisor(10, 1);
        glEnableVertexAttribArray(11);
        glVertexAttribIPointer(11, 2, GL_UNSIGNED_INT, sizeof(InstanceData), (GLvoid*)(9 * sizeof(float) + 2 * sizeof(unsigned)));
        glVertexAttribDivisor(11, 1);
        InstanceVAO[std::pair<video::E_VERTEX_TYPE, InstanceType>(tp, InstanceTypeShadow)] = vao;

        vao = createVAO(vbo[tp], ibo[tp], tp);
        glBindBuffer(GL_ARRAY_BUFFER, instance_vbo[InstanceTypeRSM]);

        glEnableVertexAttribArray(7);
        glVertexAttribPointer(7, 3, GL_FLOAT, GL_FALSE, sizeof(InstanceData), 0);
        glVertexAttribDivisor(7, 1);
        glEnableVertexAttribArray(8);
        glVertexAttribPointer(8, 3, GL_FLOAT, GL_FALSE, sizeof(InstanceData), (GLvoid*)(3 * sizeof(float)));
        glVertexAttribDivisor(8, 1);
        glEnableVertexAttribArray(9);
        glVertexAttribPointer(9, 3, GL_FLOAT, GL_FALSE, sizeof(InstanceData), (GLvoid*)(6 * sizeof(float)));
        glVertexAttribDivisor(9, 1);
        glEnableVertexAttribArray(10);
        glVertexAttribIPointer(10, 2, GL_UNSIGNED_INT, sizeof(InstanceData), (GLvoid*)(9 * sizeof(float)));
        glVertexAttribDivisor(10, 1);
        glEnableVertexAttribArray(11);
        glVertexAttribIPointer(11, 2, GL_UNSIGNED_INT, sizeof(InstanceData), (GLvoid*)(9 * sizeof(float) + 2 * sizeof(unsigned)));
        glVertexAttribDivisor(11, 1);
        InstanceVAO[std::pair<video::E_VERTEX_TYPE, InstanceType>(tp, InstanceTypeRSM)] = vao;

        vao = createVAO(vbo[tp], ibo[tp], tp);
        glBindBuffer(GL_ARRAY_BUFFER, instance_vbo[InstanceTypeGlow]);

        glEnableVertexAttribArray(7);
        glVertexAttribPointer(7, 3, GL_FLOAT, GL_FALSE, sizeof(GlowInstanceData), 0);
        glVertexAttribDivisor(7, 1);
        glEnableVertexAttribArray(8);
        glVertexAttribPointer(8, 3, GL_FLOAT, GL_FALSE, sizeof(GlowInstanceData), (GLvoid*)(3 * sizeof(float)));
        glVertexAttribDivisor(8, 1);
        glEnableVertexAttribArray(9);
        glVertexAttribPointer(9, 3, GL_FLOAT, GL_FALSE, sizeof(GlowInstanceData), (GLvoid*)(6 * sizeof(float)));
        glVertexAttribDivisor(9, 1);
        glEnableVertexAttribArray(12);
        glVertexAttribPointer(12, 4, GL_UNSIGNED_BYTE, GL_TRUE, sizeof(GlowInstanceData), (GLvoid*)(9 * sizeof(float)));
        glVertexAttribDivisor(12, 1);
        InstanceVAO[std::pair<video::E_VERTEX_TYPE, InstanceType>(tp, InstanceTypeGlow)] = vao;
        glBindVertexArray(0);
    }



}

size_t VAOManager::getVertexPitch(enum VTXTYPE tp) const
{
    switch (tp)
    {
    case VTXTYPE_STANDARD:
        return getVertexPitchFromType(video::EVT_STANDARD);
    case VTXTYPE_TCOORD:
        return getVertexPitchFromType(video::EVT_2TCOORDS);
    case VTXTYPE_TANGENT:
        return getVertexPitchFromType(video::EVT_TANGENTS);
    default:
        assert(0 && "Wrong vtxtype");
        return -1;
    }
}

VAOManager::VTXTYPE VAOManager::getVTXTYPE(video::E_VERTEX_TYPE type)
{
    switch (type)
    {
    default:
        assert(0 && "Wrong vtxtype");
    case video::EVT_STANDARD:
        return VTXTYPE_STANDARD;
    case video::EVT_2TCOORDS:
        return VTXTYPE_TCOORD;
    case video::EVT_TANGENTS:
        return VTXTYPE_TANGENT;
    }
};

void VAOManager::append(scene::IMeshBuffer *mb, VTXTYPE tp)
{
    size_t old_vtx_cnt = last_vertex[tp];
    size_t old_idx_cnt = last_index[tp];

    regenerateBuffer(tp, old_vtx_cnt + mb->getVertexCount(), old_idx_cnt + mb->getIndexCount());
    if (irr_driver->hasBufferStorageExtension())
    {
        void *tmp = (char*)VBOPtr[tp] + old_vtx_cnt * getVertexPitch(tp);
        memcpy(tmp, mb->getVertices(), mb->getVertexCount() * getVertexPitch(tp));
    }
    else
    {
        glBindBuffer(GL_ARRAY_BUFFER, vbo[tp]);
        glBufferSubData(GL_ARRAY_BUFFER, old_vtx_cnt * getVertexPitch(tp), mb->getVertexCount() * getVertexPitch(tp), mb->getVertices());
    }
    if (irr_driver->hasBufferStorageExtension())
    {
        void *tmp = (char*)IBOPtr[tp] + old_idx_cnt * sizeof(u16);
        memcpy(tmp, mb->getIndices(), mb->getIndexCount() * sizeof(u16));
    }
    else
    {
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo[tp]);
        glBufferSubData(GL_ELEMENT_ARRAY_BUFFER, old_idx_cnt * sizeof(u16), mb->getIndexCount() * sizeof(u16), mb->getIndices());
    }

    mappedBaseVertex[tp][mb] = old_vtx_cnt;
    mappedBaseIndex[tp][mb] = old_idx_cnt * sizeof(u16);
}

std::pair<unsigned, unsigned> VAOManager::getBase(scene::IMeshBuffer *mb)
{
    VTXTYPE tp = getVTXTYPE(mb->getVertexType());
    if (mappedBaseVertex[tp].find(mb) == mappedBaseVertex[tp].end())
    {
        assert(mappedBaseIndex[tp].find(mb) == mappedBaseIndex[tp].end());
        append(mb, tp);
        regenerateVAO(tp);
        regenerateInstancedVAO();
    }

    std::unordered_map<scene::IMeshBuffer*, unsigned>::iterator It;
    It = mappedBaseVertex[tp].find(mb);
    assert(It != mappedBaseVertex[tp].end());
    unsigned vtx = It->second;
    It = mappedBaseIndex[tp].find(mb);
    assert(It != mappedBaseIndex[tp].end());
    return std::pair<unsigned, unsigned>(vtx, It->second);
}
