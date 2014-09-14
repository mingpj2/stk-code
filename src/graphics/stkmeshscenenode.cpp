#include "stkmeshscenenode.hpp"
#include "stkmesh.hpp"
#include "graphics/irr_driver.hpp"
#include "graphics/glwrap.hpp"
#include "tracks/track.hpp"
#include <ISceneManager.h>
#include <IMaterialRenderer.h>
#include "config/user_config.hpp"
#include "graphics/callbacks.hpp"
#include "graphics/camera.hpp"
#include "modes/world.hpp"
#include "utils/helpers.hpp"
#include "utils/tuple.hpp"
#include "utils/cpp2011.hpp"

STKMeshSceneNode::STKMeshSceneNode(irr::scene::IMesh* mesh, ISceneNode* parent, irr::scene::ISceneManager* mgr, irr::s32 id,
    const irr::core::vector3df& position,
    const irr::core::vector3df& rotation,
    const irr::core::vector3df& scale, bool createGLMeshes) :
    CMeshSceneNode(mesh, parent, mgr, id, position, rotation, scale)
{
    isDisplacement = false;
    immediate_draw = false;
    update_each_frame = false;
    isGlow = false;

    if (createGLMeshes)
        this->createGLMeshes();
}

void STKMeshSceneNode::setReloadEachFrame(bool val)
{
    update_each_frame = val;
    if (val)
        immediate_draw = true;
}

void STKMeshSceneNode::createGLMeshes()
{
    for (u32 i = 0; i<Mesh->getMeshBufferCount(); ++i)
    {
        scene::IMeshBuffer* mb = Mesh->getMeshBuffer(i);
        GLmeshes.push_back(allocateMeshBuffer(mb));
    }
    isMaterialInitialized = false;
    isGLInitialized = false;
}

void STKMeshSceneNode::cleanGLMeshes()
{
    for (u32 i = 0; i < GLmeshes.size(); ++i)
    {
        GLMesh mesh = GLmeshes[i];
        if (!mesh.vertex_buffer)
            continue;
        if (mesh.vao)
            glDeleteVertexArrays(1, &(mesh.vao));
        if (mesh.vertex_buffer)
            glDeleteBuffers(1, &(mesh.vertex_buffer));
        if (mesh.index_buffer)
            glDeleteBuffers(1, &(mesh.index_buffer));
    }
    GLmeshes.clear();
    for (unsigned i = 0; i < MAT_COUNT; i++)
        MeshSolidMaterial[i].clearWithoutDeleting();
}

void STKMeshSceneNode::setMesh(irr::scene::IMesh* mesh)
{
    CMeshSceneNode::setMesh(mesh);
    cleanGLMeshes();
    createGLMeshes();
}

STKMeshSceneNode::~STKMeshSceneNode()
{
    cleanGLMeshes();
}

void STKMeshSceneNode::drawGlow(const GLMesh &mesh)
{
    assert(mesh.VAOType == video::EVT_STANDARD);

    GLenum ptype = mesh.PrimitiveType;
    GLenum itype = mesh.IndexType;
    size_t count = mesh.IndexCount;
    MeshShader::ColorizeShader::getInstance()->setUniforms(AbsoluteTransformation, video::SColorf(glowcolor.getRed() / 255.f, glowcolor.getGreen() / 255.f, glowcolor.getBlue() / 255.f));
    glDrawElementsBaseVertex(ptype, count, itype, (GLvoid *)mesh.vaoOffset, mesh.vaoBaseVertex);
}

void STKMeshSceneNode::updatevbo()
{
    for (unsigned i = 0; i < Mesh->getMeshBufferCount(); ++i)
    {
        scene::IMeshBuffer* mb = Mesh->getMeshBuffer(i);
        if (!mb)
            continue;
        GLMesh &mesh = GLmeshes[i];
        glDeleteBuffers(1, &(mesh.vertex_buffer));
        glDeleteBuffers(1, &(mesh.index_buffer));
        glDeleteVertexArrays(1, &(mesh.vao));

        fillLocalBuffer(mesh, mb);
        mesh.vao = createVAO(mesh.vertex_buffer, mesh.index_buffer, mb->getVertexType());
    }
}

void STKMeshSceneNode::updateNoGL()
{
    Box = Mesh->getBoundingBox();

    if (!isMaterialInitialized)
    {
        irr::video::IVideoDriver* driver = irr_driver->getVideoDriver();
        for (u32 i = 0; i < Mesh->getMeshBufferCount(); ++i)
        {
            scene::IMeshBuffer* mb = Mesh->getMeshBuffer(i);
            if (!mb)
                continue;

            video::E_MATERIAL_TYPE type = mb->getMaterial().MaterialType;
            f32 MaterialTypeParam = mb->getMaterial().MaterialTypeParam;
            video::IMaterialRenderer* rnd = driver->getMaterialRenderer(type);
            if (!isObject(type))
            {
#ifdef DEBUG
                Log::warn("material", "Unhandled (static) material type : %d", type);
#endif
                continue;
            }

            GLMesh &mesh = GLmeshes[i];
            if (rnd->isTransparent())
            {
                TransparentMaterial TranspMat = MaterialTypeToTransparentMaterial(type, MaterialTypeParam);
                if (!immediate_draw)
                    TransparentMesh[TranspMat].push_back(&mesh);
            }
            else
            {
                assert(!isDisplacement);
                MeshMaterial MatType = MaterialTypeToMeshMaterial(type, mb->getVertexType());
                if (!immediate_draw)
                    MeshSolidMaterial[MatType].push_back(&mesh);
            }
        }
        isMaterialInitialized = true;
    }

    for (u32 i = 0; i < Mesh->getMeshBufferCount(); ++i)
    {
        scene::IMeshBuffer* mb = Mesh->getMeshBuffer(i);
        if (!mb)
            continue;
        GLmeshes[i].TextureMatrix = getMaterial(i).getTextureMatrix(0);
    }
}

void STKMeshSceneNode::updateGL()
{
    if (isGLInitialized)
        return;
    for (u32 i = 0; i < Mesh->getMeshBufferCount(); ++i)
    {
        scene::IMeshBuffer* mb = Mesh->getMeshBuffer(i);
        if (!mb)
            continue;
        GLMesh &mesh = GLmeshes[i];

        irr::video::IVideoDriver* driver = irr_driver->getVideoDriver();
        video::E_MATERIAL_TYPE type = mb->getMaterial().MaterialType;
        f32 MaterialTypeParam = mb->getMaterial().MaterialTypeParam;
        video::IMaterialRenderer* rnd = driver->getMaterialRenderer(type);


        if (!rnd->isTransparent())
        {
            MeshMaterial MatType = MaterialTypeToMeshMaterial(type, mb->getVertexType()); 
            if (!immediate_draw)
                InitTextures(mesh, MatType);
        }

        if (!immediate_draw && irr_driver->hasARB_base_instance())
        {
            std::pair<unsigned, unsigned> p = VAOManager::getInstance()->getBase(mb);
            mesh.vaoBaseVertex = p.first;
            mesh.vaoOffset = p.second;
        }
        else
        {
            fillLocalBuffer(mesh, mb);
            mesh.vao = createVAO(mesh.vertex_buffer, mesh.index_buffer, mb->getVertexType());
            glBindVertexArray(0);
        }
    }
    isGLInitialized = true;
}

void STKMeshSceneNode::OnRegisterSceneNode()
{
    if (isDisplacement)
        SceneManager->registerNodeForRendering(this, scene::ESNRP_TRANSPARENT);
    else
        CMeshSceneNode::OnRegisterSceneNode();
}

static video::ITexture *spareWhiteTex = 0;

void STKMeshSceneNode::render()
{
    irr::video::IVideoDriver* driver = irr_driver->getVideoDriver();

    if (!Mesh || !driver)
        return;

    ++PassCount;

    updateNoGL();
    updateGL();

    bool isTransparent;

    for (u32 i = 0; i < Mesh->getMeshBufferCount(); ++i)
    {
        scene::IMeshBuffer* mb = Mesh->getMeshBuffer(i);
        if (!mb)
            continue;

        video::E_MATERIAL_TYPE type = mb->getMaterial().MaterialType;
        video::IMaterialRenderer* rnd = driver->getMaterialRenderer(type);

        isTransparent = rnd->isTransparent();
        break;
    }

    if ((irr_driver->getPhase() == SOLID_NORMAL_AND_DEPTH_PASS) && immediate_draw && !isTransparent)
    {
        core::matrix4 invmodel;
        AbsoluteTransformation.getInverse(invmodel);

        glDisable(GL_CULL_FACE);
        if (update_each_frame)
            updatevbo();
        glUseProgram(MeshShader::ObjectPass1Shader::getInstance()->Program);
        // Only untextured
        for (unsigned i = 0; i < GLmeshes.size(); i++)
        {
            irr_driver->IncreaseObjectCount();
            GLMesh &mesh = GLmeshes[i];
            GLenum ptype = mesh.PrimitiveType;
            GLenum itype = mesh.IndexType;
            size_t count = mesh.IndexCount;

            if (!mesh.textures[0])
                mesh.textures[0] = getUnicolorTexture(video::SColor(255, 255, 255, 255));
            compressTexture(mesh.textures[0], true);
            if (UserConfigParams::m_azdo)
            {
                if (!mesh.TextureHandles[0])
                    mesh.TextureHandles[0] = glGetTextureSamplerHandleARB(getTextureGLuint(mesh.textures[0]), MeshShader::TransparentFogShader::getInstance()->SamplersId[0]);
                if (!glIsTextureHandleResidentARB(mesh.TextureHandles[0]))
                    glMakeTextureHandleResidentARB(mesh.TextureHandles[0]);
                MeshShader::ObjectPass1Shader::getInstance()->SetTextureHandles(createVector<uint64_t>(mesh.TextureHandles[0]));
            }
            else
                MeshShader::ObjectPass1Shader::getInstance()->SetTextureUnits(std::vector < GLuint > { getTextureGLuint(mesh.textures[0]) });
            MeshShader::ObjectPass1Shader::getInstance()->setUniforms(AbsoluteTransformation, invmodel);
            assert(mesh.vao);
            glBindVertexArray(mesh.vao);
            glDrawElements(ptype, count, itype, 0);
            glBindVertexArray(0);
        }
        glEnable(GL_CULL_FACE);
        return;
    }

    if (irr_driver->getPhase() == SOLID_LIT_PASS  && immediate_draw && !isTransparent)
    {
        core::matrix4 invmodel;
        AbsoluteTransformation.getInverse(invmodel);

        glDisable(GL_CULL_FACE);
        if (update_each_frame && !UserConfigParams::m_dynamic_lights)
            updatevbo();
        if (!spareWhiteTex)
            spareWhiteTex = getUnicolorTexture(video::SColor(255, 255, 255, 255));
        glUseProgram(MeshShader::ObjectPass2Shader::getInstance()->Program);
        // Only untextured
        for (unsigned i = 0; i < GLmeshes.size(); i++)
        {
            irr_driver->IncreaseObjectCount();
            GLMesh &mesh = GLmeshes[i];
            GLenum ptype = mesh.PrimitiveType;
            GLenum itype = mesh.IndexType;
            size_t count = mesh.IndexCount;

            if (UserConfigParams::m_azdo)
            {
                GLuint64 DiffuseHandle = glGetTextureSamplerHandleARB(irr_driver->getRenderTargetTexture(RTT_DIFFUSE), MeshShader::ObjectPass2Shader::getInstance()->SamplersId[0]);
                if (!glIsTextureHandleResidentARB(DiffuseHandle))
                    glMakeTextureHandleResidentARB(DiffuseHandle);

                GLuint64 SpecularHandle = glGetTextureSamplerHandleARB(irr_driver->getRenderTargetTexture(RTT_SPECULAR), MeshShader::ObjectPass2Shader::getInstance()->SamplersId[1]);
                if (!glIsTextureHandleResidentARB(SpecularHandle))
                    glMakeTextureHandleResidentARB(SpecularHandle);

                GLuint64 SSAOHandle = glGetTextureSamplerHandleARB(irr_driver->getRenderTargetTexture(RTT_HALF1_R), MeshShader::ObjectPass2Shader::getInstance()->SamplersId[2]);
                if (!glIsTextureHandleResidentARB(SSAOHandle))
                    glMakeTextureHandleResidentARB(SSAOHandle);

                if (!mesh.TextureHandles[0])
                    mesh.TextureHandles[0] = glGetTextureSamplerHandleARB(getTextureGLuint(spareWhiteTex), MeshShader::TransparentFogShader::getInstance()->SamplersId[0]);
                if (!glIsTextureHandleResidentARB(mesh.TextureHandles[0]))
                    glMakeTextureHandleResidentARB(mesh.TextureHandles[0]);
                MeshShader::ObjectPass2Shader::getInstance()->SetTextureHandles(createVector<uint64_t>(DiffuseHandle, SpecularHandle, SSAOHandle, mesh.TextureHandles[0]));
            }
            else
                MeshShader::ObjectPass2Shader::getInstance()->SetTextureUnits(createVector<GLuint>(
                irr_driver->getRenderTargetTexture(RTT_DIFFUSE),
                irr_driver->getRenderTargetTexture(RTT_SPECULAR),
                irr_driver->getRenderTargetTexture(RTT_HALF1_R),
                getTextureGLuint(spareWhiteTex)));
            MeshShader::ObjectPass2Shader::getInstance()->setUniforms(AbsoluteTransformation, mesh.TextureMatrix);
            assert(mesh.vao);
            glBindVertexArray(mesh.vao);
            glDrawElements(ptype, count, itype, 0);
            glBindVertexArray(0);
        }
        glEnable(GL_CULL_FACE);
        return;
    }

    if (irr_driver->getPhase() == GLOW_PASS)
    {
        glUseProgram(MeshShader::ColorizeShader::getInstance()->Program);
        for (u32 i = 0; i < Mesh->getMeshBufferCount(); ++i)
        {
            scene::IMeshBuffer* mb = Mesh->getMeshBuffer(i);
            if (!mb)
                continue;
            if (irr_driver->hasARB_base_instance())
                glBindVertexArray(VAOManager::getInstance()->getVAO(video::EVT_STANDARD));
            else
                glBindVertexArray(GLmeshes[i].vao);
            drawGlow(GLmeshes[i]);
        }
    }

    if (irr_driver->getPhase() == TRANSPARENT_PASS && isTransparent)
    {
        ModelViewProjectionMatrix = computeMVP(AbsoluteTransformation);

        if (immediate_draw)
        {
            if (update_each_frame)
                updatevbo();
            glBlendFunc(GL_ONE, GL_ONE_MINUS_SRC_ALPHA);

            if (World::getWorld() && World::getWorld()->isFogEnabled())
            {
                glUseProgram(MeshShader::TransparentFogShader::getInstance()->Program);
                for (unsigned i = 0; i < GLmeshes.size(); i++)
                {
                    GLMesh &mesh = GLmeshes[i];
                    irr_driver->IncreaseObjectCount();
                    GLenum ptype = mesh.PrimitiveType;
                    GLenum itype = mesh.IndexType;
                    size_t count = mesh.IndexCount;

                    const Track * const track = World::getWorld()->getTrack();

                    // This function is only called once per frame - thus no need for setters.
                    const float fogmax = track->getFogMax();
                    const float startH = track->getFogStartHeight();
                    const float endH = track->getFogEndHeight();
                    const float start = track->getFogStart();
                    const float end = track->getFogEnd();
                    const video::SColor tmpcol = track->getFogColor();

                    video::SColorf col(tmpcol.getRed() / 255.0f,
                        tmpcol.getGreen() / 255.0f,
                        tmpcol.getBlue() / 255.0f);

                    if (!mesh.textures[0])
                        mesh.textures[0] = getUnicolorTexture(video::SColor(255, 255, 255, 255));
                    compressTexture(mesh.textures[0], true);
                    if (UserConfigParams::m_azdo)
                    {
                        if (!mesh.TextureHandles[0])
                            mesh.TextureHandles[0] = glGetTextureSamplerHandleARB(getTextureGLuint(mesh.textures[0]), MeshShader::TransparentFogShader::getInstance()->SamplersId[0]);
                        if (!glIsTextureHandleResidentARB(mesh.TextureHandles[0]))
                            glMakeTextureHandleResidentARB(mesh.TextureHandles[0]);
                        MeshShader::TransparentFogShader::getInstance()->SetTextureHandles(createVector<uint64_t>(mesh.TextureHandles[0]));
                    }
                    else
                        MeshShader::TransparentFogShader::getInstance()->SetTextureUnits(std::vector<GLuint>{ getTextureGLuint(mesh.textures[0]) });
                    MeshShader::TransparentFogShader::getInstance()->setUniforms(AbsoluteTransformation, mesh.TextureMatrix, fogmax, startH, endH, start, end, col);

                    assert(mesh.vao);
                    glBindVertexArray(mesh.vao);
                    glDrawElements(ptype, count, itype, 0);
                    glBindVertexArray(0);
                }
            }
            else
            {
                glUseProgram(MeshShader::TransparentShader::getInstance()->Program);
                for (unsigned i = 0; i < GLmeshes.size(); i++)
                {
                    irr_driver->IncreaseObjectCount();
                    GLMesh &mesh = GLmeshes[i];
                    GLenum ptype = mesh.PrimitiveType;
                    GLenum itype = mesh.IndexType;
                    size_t count = mesh.IndexCount;
                    if (!mesh.textures[0])
                        mesh.textures[0] = getUnicolorTexture(video::SColor(255, 255, 255, 255));
                    compressTexture(mesh.textures[0], true);
                    if (UserConfigParams::m_azdo)
                    {
                        if (!mesh.TextureHandles[0])
                            mesh.TextureHandles[0] = glGetTextureSamplerHandleARB(getTextureGLuint(mesh.textures[0]), MeshShader::TransparentShader::getInstance()->SamplersId[0]);
                        if (!glIsTextureHandleResidentARB(mesh.TextureHandles[0]))
                            glMakeTextureHandleResidentARB(mesh.TextureHandles[0]);
                        MeshShader::TransparentShader::getInstance()->SetTextureHandles(createVector<uint64_t>(mesh.TextureHandles[0]));
                    }
                    else
                        MeshShader::TransparentShader::getInstance()->SetTextureUnits(std::vector<GLuint>{ getTextureGLuint(mesh.textures[0]) });

                    MeshShader::TransparentShader::getInstance()->setUniforms(AbsoluteTransformation, mesh.TextureMatrix);
                    assert(mesh.vao);
                    glBindVertexArray(mesh.vao);
                    glDrawElements(ptype, count, itype, 0);
                    glBindVertexArray(0);
                }
            }
            return;
        }

        GLMesh* mesh;
        if (!TransparentMesh[TM_BUBBLE].empty())
            glUseProgram(MeshShader::BubbleShader::Program);
        if (irr_driver->hasARB_base_instance())
            glBindVertexArray(VAOManager::getInstance()->getVAO(video::EVT_STANDARD));
        for_in(mesh, TransparentMesh[TM_BUBBLE])
        {
            if (irr_driver->hasARB_base_instance())
                glBindVertexArray(mesh->vao);
            drawBubble(*mesh, ModelViewProjectionMatrix);
        }
        return;
    }
}
