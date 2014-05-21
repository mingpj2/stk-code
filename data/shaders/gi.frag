//----------------------------------------------------------------------------------
//
//  Radiance hints GI reconstruction fragment shader.
//  
//  G. Papaioannou (gepap@aueb.gr), 2011
//
//
//----------------------------------------------------------------------------------
//
// Abreviations: 
// CSS: Canonical Screen Space
// WCS: World Coordinate System
// ECS: Eye (camera) Coordinate System
//
//----------------------------------------------------------------------------------
//
// This shader assumes that a GI buffer-covering quad is being drawn. The vertex shader
// only passes the tex. coords through to the fragment shader. Tex coord set 0 is used
// for recovering the CSS coordinates of the current fragment. Bottom left: (0,0), top
// right: (1,1).  
//
// Notes:
// This is the 2nd order SH (4 coef) version of the demo/example shader source code. 
// higher order SHs are similarly utilized. The shader could also be further optimized 
// in order to avoid re-computing certain quantities.


// ----------------  INPUT ---------------------------------------------------------

// Deferred renderer normal and depth buffers:
uniform sampler2D ntex;
uniform sampler2D dtex;

uniform sampler3D SHR;
uniform sampler3D SHG;
uniform sampler3D SHB;

uniform float R_wcs = 10.;          // Rmax: maximum sampling distance (in WCS units)
uniform float factor = 1.;         // GI contribution multiplier 
uniform vec3 extents;          // Bounding box limits of the radiance hints volume
uniform sampler3D points[4];  // Radiance Hint buffers:
                              // points[0]: dist_min, dist_max, dist_ave, 1.0
                              // points[1]: Lr(1,1) Lr(1,-1) Lr(1,0) Lr(0,0) 
                              // points[2]: Lg(1,1) Lg(1,-1) Lg(1,0) Lg(0,0) 
                              // points[3]: Lb(1,1) Lb(1,-1) Lb(1,0) Lb(0,0) 

layout (std140) uniform MatrixesData
{
    mat4 ViewMatrix;
    mat4 ProjectionMatrix;
    mat4 InverseViewMatrix;
    mat4 InverseProjectionMatrix;
    mat4 ShadowViewProjMatrixes[4];
    vec2 screen;
};

// ----------------  SH functions -------------------------------------------------

vec4 SHBasis (const in vec3 dir)
{ 
    float   L00  = 0.282095;
    float   L1_1 = 0.488603 * dir.y;
    float   L10  = 0.488603 * dir.z;
    float   L11  = 0.488603 * dir.x;
    return vec4 (L11, L1_1, L10, L00);
}

vec3 SH2RGB (in vec4 sh_r, in vec4 sh_g, in vec4 sh_b, in vec3 dir)
{
    vec4 Y = vec4(1.023326*dir.x, 1.023326*dir.y, 1.023326*dir.z, 0.886226);
    return vec3 (dot(Y,sh_r), dot(Y,sh_g), dot(Y,sh_b));
}

// ----------------  Coordinate transformations -----------------------------------

vec3 VectorECS2WCS(in vec3 pos)
{
    vec4 vector_WCS = transpose(ViewMatrix) * vec4(pos,0);
    return vector_WCS.xyz;
}

// ----------------  Main shader function -----------------------------------------

in vec2 uv;
out vec4 Diffuse;
out vec4 Specular;

vec3 DecodeNormal(vec2 n);
vec4 getPosFromUVDepth(vec3 uvDepth, mat4 InverseProjectionMatrix);

void main(void)
{
    // Accumulated global illumination, initialized at (0,0,0)
    vec3 GI = vec3(0.);

    // Blending factor for current GI estimation frame
    // If blending < 1, current result is blended with previous GI values (temporal smoothing)
    float blending = 0.8;

    float depth = texture2D(dtex, uv).x;
    // Discard background fragments
    if (depth==1.0) discard;

    ivec3 sz = textureSize(points[0],0);

    // convert fragment position and normal to WCS
    vec3 pos_wcs = (InverseViewMatrix * getPosFromUVDepth(vec3(uv, depth), InverseProjectionMatrix)).xyz;
    vec3 normal_ecs = normalize(DecodeNormal(2. * texture(ntex, uv).xy - 1.));
    vec3 normal_wcs = normalize(VectorECS2WCS(normal_ecs));

    // determine volume texture coordinate of current fragment location
    vec3 uvw = pos_wcs / extents;

    float denom = 0.05;

    // Sample the RH volume at 4 locations, one directly above the shaded point,
    // three on a ring 80degs away from the normal direction. All samples are
    // moved away from the shaded point to avoid sampling RHs behind the surface.
    // You can introduce a random rotation of the samples around the normal
    vec3 rnd = vec3(0,0,0);

    // Generate the sample locations
    vec3 v_rand = vec3(0.5);
    vec3 v_1 = normalize(cross(normal_wcs,v_rand));
    vec3 v_2 = cross(normal_wcs,v_1);
    vec3 D[4];
    D[0] = vec3(1.0,0.0,0.0);
    int i;
    for (i=1; i<4; i++)
    {
        D[i] = vec3(0.1, 0.8*cos((rnd.x*1.5+i)*6.2832/3.0), 0.8*sin((rnd.x*1.5+i)*6.2832/3.0));
        D[i] = normalize(D[i]);
    }

    for (i=0; i<4; i++)
    {
        vec3 sdir = normal_wcs*D[i].x + v_1*D[i].y + v_2*D[i].z;
        vec3 uvw_new = 0.5*normal_wcs/sz+ sdir/sz + uvw;

        vec4 rh_shr    = texture3D(points[1],uvw_new);
        vec4 rh_shg    = texture3D(points[2],uvw_new);
        vec4 rh_shb    = texture3D(points[3],uvw_new);

        vec3 rh_pos    = extents * uvw_new;
        vec3 path = rh_pos - pos_wcs;
        float dist = length(path);
        float rel_dist = dist/R_wcs;
        dist/=length(extents/sz);
        path = normalize(path);
        float contrib = dist>0.005?1.0:0.0;
        GI+= contrib * SH2RGB (rh_shr, rh_shg, rh_shb, -normal_wcs);
        denom+=contrib;
    }
    GI *= factor / denom;

    // Note: tone mapping is normally required on the final GI buffer
    Diffuse = vec4(GI,blending);
    Specular = vec4(0.);
}
