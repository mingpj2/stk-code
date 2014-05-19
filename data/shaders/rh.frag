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

// ----------------  FLAGS ---------------------------------------------------------

// #define DEPTH_OCCL    // if defined, depth-based RSM sample occlusion is enabled.

// ----------------  INPUT ---------------------------------------------------------

uniform ivec3 resolution;       // The volume resolution 
uniform int slice;              // The current volume slice
uniform sampler3D Noise;        // A pre-computed 3D noise texture (32X32X32). Value range (r,g,b): [0,1]
uniform float R_wcs;            // Rmax: maximum sampling distance (in WCS units)
uniform vec3 bbox_min;          // Bounding box limits of the radiance hints volume
uniform vec3 bbox_max;          // 
uniform sampler2D Shadow;       // RSM depth
uniform sampler2D ShadowColor;  // RSM vpl flux
uniform sampler2D ShadowNormal; // RSM normals
uniform sampler2D Depth;        // camera depth buffer
uniform int samples;            // Number of RSM samples
uniform mat4 L;                 // Light transformation matrix ( WCS -> LCS(RSM))
uniform mat4 L_inv;             // Inverse light transformation matrix ( LCS(RSM) -> WCS)
uniform mat4 L_ecs;             // Light modelview transformation matrix ( WCS -> Light ECS )
uniform mat4  MVP;              // Final modelview and projection camera matrix (CSS -> WCS)    
uniform int num_lights;         // Number of GI lights
uniform float spread;           // RSM parametric sampling radius
uniform vec3 light_pos;         // Current light position in WCS coords
uniform vec3 light_dir;         // Current light direction in WCS coords


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

vec3 PointWCS2CSS(in vec3 sample) 
{ 
    vec4 p_css = MVP*vec4(sample,1); 
    return p_css.xyz/p_css.w; 
} 

vec3 VectorLECS2WCS(in vec3 sample) 
{ 
    // For vectors, the transposed inverse matrix from Light ECS to WCS is used 
    // (i.e. the transposed L_ecs). You could also pass the matrix transposed 
    // outside the shader.
    return (transpose(L_ecs)*vec4(sample,0)).xyz;
} 

vec2 ShadowProjection( in vec3 point_WCS ) 
{ 
    // note: projected points on the RSM are clamped to the map extents,
    // not rejected
    vec4 pos_LCS = L*vec4(point_WCS+vec3(0.01,0.01,0.01),1.0); 
    pos_LCS/=pos_LCS.w; 
    float reverse = sign(dot(light_dir,point_WCS-light_pos)); 
    vec2 uv = vec2( reverse*0.5*pos_LCS.x + 0.5, reverse*0.5*pos_LCS.y + 0.5); 
    return clamp(uv,vec2(0.0,0.0),vec2(1.0,1.0)); 
} 


// ----------------  Main shader function -----------------------------------------


void main(void) 
{ 
    // Determine the RH position (WCS) 
    int   gy = gl_FragCoord.y; 
    int   gx = gl_FragCoord.x; 
    int   gz = slice; 
    vec3  extents = bbox_max-bbox_min; 
    vec3  stratum; 
    stratum.x = extents.x/(resolution.x-1); 
    stratum.y = extents.y/(resolution.y-1); 
    stratum.z = extents.z/(resolution.z-1); 
    vec3  pos = bbox_min + vec3(gx, gy, gz) * stratum; 

    // Project the RH location on the RSM and determine the 
    // center of the sampling disk
    vec2 l_uv = ShadowProjection(pos); 
  
  
    // initialize measured distances from RH location to the RSM samples
    vec4 SH_dist_ave = vec4(0,0,0,0); 
    float dist, dist_min=R_wcs, dist_max=0.0, dist_ave=0.0; 

    // Declare and initialize various parameters
    vec4  smc, smp, smn; 
    vec4  SHr=vec4(0.0,0.0,0.0,0.0); // accumulated SH coefs for the incident radiance from 
    vec4  SHg=vec4(0.0,0.0,0.0,0.0); // all RSM samples
    vec4  SHb=vec4(0.0,0.0,0.0,0.0); 
    vec3  normal; 
    vec4 color; 
 
    for (int i=0;i<samples;i++) 
    { 
	// produce a new sample location on the RSM texture
        vec3 rnd = 2.0*texture3D(Noise, 14*pos/extents+vec3(i,0,0)/samples).xyz-vec3(1.0,1.0,1.0); 
	vec2 uv = l_uv+vec2(rnd.x*spread*cos(6.283*rnd.y),rnd.x*spread*sin(6.283*rnd.y)); 
	
        // produce a new sampling location in the RH stratum
	vec3 p = pos+(0.5*rnd)*stratum;

	// determine the the WCS coordinates of the RSM sample position (smp) and normal (smn)
        // and read the corresponding lighting (smc)
        float depth = (texture2D(Shadow,uv)).z; 
	vec4 pos_LCS = vec4 (uv.x*2.0-1, uv.y*2.0-1.0, depth*2.0-1.0,1.0); 
	vec4 isect4 = L_inv*pos_LCS; 
	smp = isect4.xyz/isect4.w; 
	smc=texture2D(ShadowColor,uv); 
	vec4 ntex = texture2D(ShadowNormal,uv); 
	smn = vec3(ntex.x*2.0-1.0,ntex.y*2.0-1.0,ntex.z); 
	smn = normalize(VectorLECS2WCS(normalize(smn))); 
	
	// Normalize distance to RSM sample
        dist = distance(p,smp)/R_wcs; 
        // Determine the incident direction. 
        // Avoid very close samples (and numerical instability problems)
        vec3 dir = dist<=0.007?vec3(0,0,0):normalize(p-smp);
	float dotprod = max(dot(dir,smn),0.0); 
	FF = dotprod/(0.1+dist*dist);
	
#ifdef DEPTH_OCCL 
        // ---- Depth-buffer-based RSM sample occlusion

	// Initialize visibility to 1
        float depth_visibility = 1.0; 

	// set number of visibility samples along the line of sight. Can be set with #define
        float vis_samples = 4.0; // 3 to 8
	vec3 Qj; 
	vec3 Qcss; 
	
        // Check if the current RH point is hidden from view. If it is, then "visible" line-of-sight 
        // samples should also be hidden from view. If the RH point is visible to the camera, the 
        // same should hold for the visibility samples in order not to attenuate the light.
        Qcss = PointWCS2CSS(p); 
	float rh_visibility = Qcss.z<(2.0*texture2D(Depth,0.5*Qcss.xy+vec2(0.5,0.5)).r-1.0)*1.1?1.0:-1.0; 

	// Estimate attenuation along line of sight
        for (int j=1; j<vis_samples; j++) 
	{ 
	    // determine next point along the line of sight
            Qj = smp+(j/vis_samples)*(pos-smp); 
	    Qcss = PointWCS2CSS(Qj); 
            // modulate the visibility according to the number of hidden LoS samples
	    depth_visibility -= rh_visibility*Qcss.z<rh_visibility*(2.0*texture2D(Depth,0.5*Qcss.xy+vec2(0.5,0.5)).r-1.0)?0.0:1.0/vis_samples; 
	} 
	depth_visibility = clamp(depth_visibility,0,1);
	FF *= depth_visibility;

#endif 

	color = smc*FF;
	vec4 shr, shg, shb; 
        // encode radiance into SH coefs and accumulate to RH
	RGB2SH (dir, color.rgb, shr, shg, shb); 
	SHr+=shr; 
	SHg+=shg; 
	SHb+=shb; 
        // update distance measurements
	dist_max=(dist>dist_max )?dist:dist_max; 
	dist_min=(dist<dist_min )?dist:dist_min; 
	dist_ave+=dist; 
} 

// cast samples to float to resolve some weird compiler issue. 
// For some reason, both floatXint and intXfloat are auto-boxed to int!
SHr/=(3.14159*float(samples)); 
SHg/=float(3.14159*float(samples)); 
SHb/=float(3.14159*float(samples)); 
dist_ave/=float(samples); 

gl_FragData[0] = vec4 (dist_min,dist_max,dist_ave,1.0)/num_lights;
// distances must be divided by the number of lights this shader is run for, because
// because they should not be accumulated. Ideally, the MIN/MAX frame buffer operator should
// be used instead of the ADD operator for this data channel. In this case (e.g. MIN), 
// frag data[0] should contain: dist_min, r_max-dist_max, dist_ave (not used actually), 1.0 )

gl_FragData[1] = SHr; 
gl_FragData[2] = SHg; 
gl_FragData[3] = SHb; 

}


