//----------------------------------------------------------------------------------
//
//  Radiance hints generation fragment shader (first bounce: RSM sampling).
//  
//  G. Papaioannou (gepap@aueb.gr), 2011
//
//----------------------------------------------------------------------------------
//
// Abreviations: 
// CSS: Canonical Screen Space
// WCS: World Coordinate System
// ECS: Eye (camera) Coordinate System
// LCS: Light Clip Space Coordinate System
//
//----------------------------------------------------------------------------------
//
// The shader is executed for every RH (voxel) of the RH buffer (3D texture), i.e.
// for every fragment of a sweep plane through the volume. For every slice, a 
// quad is drawn. The process is repeated for every RSM that contributes to the GI. 
// For this reason, the blending equation should be set to ADD and blending coefs to
// GL_ONE (src), GL_ONE (tgt)
// If the blending parameters for color attachment 0 (distance data) can be separately 
// set, use MIN as the blending equation.
//
// Notes:
// This is the 2nd order SH (4 coef) version of the demo/example shader source code. 
// higher order SHs are similarly utilized. The shader could also be further optimized 
// in order to avoid re-computing certain quantities.

// ----------------  INPUT ---------------------------------------------------------

uniform int slice;              // The current volume slice
uniform float R_wcs;            // Rmax: maximum sampling distance (in WCS units)
uniform vec3 extents;          // Bounding box limits of the radiance hints volume
uniform sampler2D dtex;         // RSM depth
uniform sampler2D ctex;         // RSM vpl flux
uniform sampler2D ntex;         // RSM normals
uniform int num_lights;         // Number of GI lights

layout (std140) uniform MatrixesData
{
    mat4 ViewMatrix;
    mat4 ProjectionMatrix;
    mat4 InverseViewMatrix;
    mat4 InverseProjectionMatrix;
    mat4 ShadowViewProjMatrixes[4];
    vec2 screen;
};

//out vec4 VariousData;
out vec4 SHRed;
out vec4 SHGreen;
out vec4 SHBlue;

#define SAMPLES 6
#define DIM_X 32
#define DIM_Y 32
#define DIM_Z 32

// ----------------  SH functions -------------------------------------------------

vec4 SHBasis (const in vec3 dir)
{
    float   L00  = 0.282095;
    float   L1_1 = 0.488603 * dir.y;
    float   L10  = 0.488603 * dir.z;
    float   L11  = 0.488603 * dir.x;
    return vec4 (L11, L1_1, L10, L00);
}

void RGB2SH (in vec3 dir, in vec3 L, out vec4 sh_r, out vec4 sh_g, out vec4 sh_b)
{
    vec4 sh = SHBasis (dir);
    sh_r = L.r * sh;
    sh_g = L.g * sh;
    sh_b = L.b * sh;
}

// ----------------  Coordinate transformations -----------------------------------

vec3 VectorLECS2WCS(vec3 pos)
{
    // ntex store normals in world space
    return pos;
}

vec2 ShadowProjection(vec3 pos)
{
    vec4 shadowcoord = ShadowViewProjMatrixes[2] * vec4(pos, 1.0);
    shadowcoord /= shadowcoord.w;
    return shadowcoord.xy * 0.5 + 0.5;
}


// ----------------  Main shader function -----------------------------------------

const vec3 rand_samples[6] = {
    vec3(1., 0., 0.),
    vec3(-1., 0., 0.),
    vec3(0., 1., 0.),
    vec3(0., -1., 0.),
    vec3(0., 0., 1.),
    vec3(0., 0., -1.),
};

void main(void)
{
    // Determine the RH center
    int   gy = int(gl_FragCoord.y);
    int   gx = int(gl_FragCoord.x);
    int   gz = slice;
    vec3  stratum = extents / vec3(DIM_X - 1., DIM_Y - 1., DIM_Z - 1.);
    vec3  RHcenter = vec3(gx, gy, gz) * stratum;

    // Project the RH location on the RSM and determine the
    // center of the sampling disk
    vec2 l_uv = ShadowProjection(RHcenter);

    // initialize measured distances from RH location to the RSM samples
    vec4 SH_dist_ave = vec4(0,0,0,0);
    float dist, dist_min=R_wcs, dist_max=0.0, dist_ave=0.0;

    // Declare and initialize various parameters
    vec3  smc, smp, smn;
    vec4  SHr = vec4(0.); // accumulated SH coefs for the incident radiance from
    vec4  SHg = vec4(0.); // all RSM samples
    vec4  SHb = vec4(0.);
    vec3 color;

    for (int i = 0; i < SAMPLES; i++)
    {
        // produce a new sample location on the RSM texture
        vec3 rnd = rand_samples[i] / 2.;
        vec2 uv = l_uv + vec2(rnd.x * cos(6.283*rnd.y), rnd.x * sin(6.283 * rnd.y));
        uv /= screen;

        // produce a new sampling location in the RH stratum
        vec3 p = RHcenter + (0.5 * rnd) * stratum;

        // determine the the WCS coordinates of the RSM sample position (smp) and normal (smn)
        // and read the corresponding lighting (smc)
        float depth = texture2D(dtex, uv).z;
        vec4 pos_LCS = inverse(ShadowViewProjMatrixes[2]) * (2. * vec4(uv, depth, 1.) - 1.);
        smp = pos_LCS.xyz / pos_LCS.w;
        smc = texture(ctex, uv).xyz;
        vec4 normal = texture(ntex, uv);
        smn = vec3(2. * normal.xy - 1., normal.z);
        smn = normalize(VectorLECS2WCS(normalize(smn)));

        // Normalize distance to RSM sample
        dist = distance(p,smp)/R_wcs;
        // Determine the incident direction.
        // Avoid very close samples (and numerical instability problems)
        vec3 dir = (dist <= 0.007) ? vec3(0.) : normalize(p-smp);
        float dotprod = max(dot(dir, smn), 0.);
        float FF = dotprod / (0.1 + dist * dist);

        color = smc * FF;
        vec4 shr, shg, shb;
        // encode radiance into SH coefs and accumulate to RH
        RGB2SH (dir, color.rgb, shr, shg, shb);
        SHr += shr;
        SHg += shg;
        SHb += shb;
        // update distance measurements
        dist_max = (dist>dist_max ) ? dist : dist_max;
        dist_min = (dist<dist_min ) ? dist : dist_min;
        dist_ave += dist;
    }

    // cast samples to float to resolve some weird compiler issue.
    // For some reason, both floatXint and intXfloat are auto-boxed to int!
    SHr /= (3.14159*float(SAMPLES));
    SHg /= float(3.14159*float(SAMPLES));
    SHb /= float(3.14159*float(SAMPLES));
    dist_ave /= float(SAMPLES);

//    VariousData = vec4(dist_min, dist_max, dist_ave, 1.) / num_lights;
    // distances must be divided by the number of lights this shader is run for, because
    // because they should not be accumulated. Ideally, the MIN/MAX frame buffer operator should
    // be used instead of the ADD operator for this data channel. In this case (e.g. MIN),
    // frag data[0] should contain: dist_min, r_max-dist_max, dist_ave (not used actually), 1.0 )

    SHRed = vec4(0.);//SHr;
    SHGreen = SHg;
    SHBlue = SHb;
}