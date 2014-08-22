#ifndef GLWRAP_HEADER_H
#define GLWRAP_HEADER_H

#include "gl_headers.hpp"

#include <vector>
#include "irr_driver.hpp"
#include "utils/log.hpp"

// already includes glext.h, which defines useful GL constants.
// COpenGLDriver has already loaded the extension GL functions we use (e.g glBeginQuery)
#include "../../lib/irrlicht/source/Irrlicht/COpenGLDriver.h"


void initGL();
GLuint LoadTFBProgram(const char * vertex_file_path, const char **varyings, unsigned varyingscount);
video::ITexture* getUnicolorTexture(video::SColor c);
void setTexture(unsigned TextureUnit, GLuint TextureId, GLenum MagFilter, GLenum MinFilter, bool allowAF = false);
GLuint LoadShader(const char * file, unsigned type);

template<typename ... Types>
void loadAndAttach(GLint ProgramID)
{
    return;
}

template<typename ... Types>
void loadAndAttach(GLint ProgramID, GLint ShaderType, const char *filepath, Types ... args)
{
    GLint ShaderID = LoadShader(filepath, ShaderType);
    glAttachShader(ProgramID, ShaderID);
    glDeleteShader(ShaderID);
    loadAndAttach(ProgramID, args...);
}

template<typename ...Types>
void printFileList()
{
    return;
}

template<typename ...Types>
void printFileList(GLint ShaderType, const char *filepath, Types ... args)
{
    Log::error("GLWrapp", filepath);
    printFileList(args...);
}

template<typename ... Types>
GLint LoadProgram(Types ... args)
{
    GLint ProgramID = glCreateProgram();
    loadAndAttach(ProgramID, args...);
    if (irr_driver->getGLSLVersion() < 330)
    {
        glBindAttribLocation(ProgramID, 0, "Position");
        glBindAttribLocation(ProgramID, 1, "Normal");
        glBindAttribLocation(ProgramID, 2, "Color");
        glBindAttribLocation(ProgramID, 3, "Texcoord");
        glBindAttribLocation(ProgramID, 4, "SecondTexcoord");
        glBindAttribLocation(ProgramID, 5, "Tangent");
        glBindAttribLocation(ProgramID, 6, "Bitangent");
        glBindAttribLocation(ProgramID, 7, "Origin");
        glBindAttribLocation(ProgramID, 8, "Orientation");
        glBindAttribLocation(ProgramID, 9, "Scale");
    }
    glLinkProgram(ProgramID);

    GLint Result = GL_FALSE;
    int InfoLogLength;
    glGetProgramiv(ProgramID, GL_LINK_STATUS, &Result);
    if (Result == GL_FALSE) {
        Log::error("GLWrapp", "Error when linking these shaders :");
        printFileList(args...);
        glGetProgramiv(ProgramID, GL_INFO_LOG_LENGTH, &InfoLogLength);
        char *ErrorMessage = new char[InfoLogLength];
        glGetProgramInfoLog(ProgramID, InfoLogLength, NULL, ErrorMessage);
        Log::error("GLWrapp", ErrorMessage);
        delete[] ErrorMessage;
    }

    GLenum glErr = glGetError();
    if (glErr != GL_NO_ERROR)
    {
        Log::warn("IrrDriver", "GLWrap : OpenGL error %i\n", glErr);
    }

    return ProgramID;
}

class GPUTimer;

class ScopedGPUTimer
{
public:
    ScopedGPUTimer(GPUTimer &);
    ~ScopedGPUTimer();
};

class GPUTimer
{
    friend class ScopedGPUTimer;
    GLuint query;
    bool initialised;
public:
    GPUTimer();
    unsigned elapsedTimeus();
};

class FrameBuffer
{
private:
    GLuint fbo;
    std::vector<GLuint> RenderTargets;
    GLuint DepthTexture;
    size_t width, height;
public:
    FrameBuffer();
    FrameBuffer(const std::vector <GLuint> &RTTs, size_t w, size_t h, bool layered = false);
    FrameBuffer(const std::vector <GLuint> &RTTs, GLuint DS, size_t w, size_t h, bool layered = false);
    ~FrameBuffer();
    void Bind();
    const std::vector<GLuint> &getRTT() const { return RenderTargets; }
    GLuint &getDepthTexture() { assert(DepthTexture); return DepthTexture; }
    size_t getWidth() const { return width; }
    size_t getHeight() const { return height; }
    static void Blit(const FrameBuffer &Src, FrameBuffer &Dst, GLbitfield mask = GL_COLOR_BUFFER_BIT, GLenum filter = GL_NEAREST);
    void BlitToDefault(size_t, size_t, size_t, size_t);
};

// core::rect<s32> needs these includes
#include <rect.h>
#include "utils/vec3.hpp"

GLuint getTextureGLuint(irr::video::ITexture *tex);
GLuint getDepthTexture(irr::video::ITexture *tex);
void resetTextureTable();
void compressTexture(irr::video::ITexture *tex, bool srgb, bool premul_alpha = false);
bool loadCompressedTexture(const std::string& compressed_tex);
void saveCompressedTexture(const std::string& compressed_tex);

class VAOManager : public Singleton<VAOManager>
{
    enum VTXTYPE { VTXTYPE_STANDARD, VTXTYPE_TCOORD, VTXTYPE_TANGENT, VTXTYPE_COUNT };
    GLuint vbo[VTXTYPE_COUNT], ibo[VTXTYPE_COUNT], vao[VTXTYPE_COUNT];
    std::vector<scene::IMeshBuffer *> storedCPUBuffer[VTXTYPE_COUNT];
    void *vtx_mirror[VTXTYPE_COUNT], *idx_mirror[VTXTYPE_COUNT];
    size_t vtx_cnt[VTXTYPE_COUNT], idx_cnt[VTXTYPE_COUNT];
    std::map<scene::IMeshBuffer*, unsigned> mappedBaseVertex[VTXTYPE_COUNT], mappedBaseIndex[VTXTYPE_COUNT];

    void regenerateBuffer(enum VTXTYPE);
    void regenerateVAO(enum VTXTYPE);
    size_t getVertexPitch(enum VTXTYPE) const;
    VTXTYPE getVTXTYPE(video::E_VERTEX_TYPE type);
    void append(scene::IMeshBuffer *, VTXTYPE tp);
public:
    VAOManager();
    std::pair<unsigned, unsigned> getBase(scene::IMeshBuffer *);
    unsigned getVBO(video::E_VERTEX_TYPE type) { return vbo[getVTXTYPE(type)]; }
    unsigned getVAO(video::E_VERTEX_TYPE type) { return vao[getVTXTYPE(type)]; }
    ~VAOManager();
};

void draw3DLine(const core::vector3df& start,
    const core::vector3df& end, irr::video::SColor color);

void draw2DImageFromRTT(GLuint texture, size_t texture_w, size_t texture_h,
    const core::rect<s32>& destRect,
    const core::rect<s32>& sourceRect, const core::rect<s32>* clipRect,
    const video::SColor &colors, bool useAlphaChannelOfTexture);

void draw2DImage(const irr::video::ITexture* texture, const irr::core::rect<s32>& destRect,
    const irr::core::rect<s32>& sourceRect, const irr::core::rect<s32>* clipRect,
    const irr::video::SColor &color, bool useAlphaChannelOfTexture);

void draw2DImage(const irr::video::ITexture* texture, const irr::core::rect<s32>& destRect,
    const irr::core::rect<s32>& sourceRect, const irr::core::rect<s32>* clipRect,
    const irr::video::SColor* const colors, bool useAlphaChannelOfTexture);

void GL32_draw2DRectangle(irr::video::SColor color, const irr::core::rect<s32>& position,
    const irr::core::rect<s32>* clip = 0);

bool hasGLExtension(const char* extension);

#endif
