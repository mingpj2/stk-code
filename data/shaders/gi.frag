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

uniform sampler2D RT_normals; // Deferred renderer normal and depth buffers:
uniform sampler2D RT_depth;   //

uniform sampler3D Noise;      // A pre-computed 3D noise texture (32X32X32). Value range (r,g,b): [0,1]
uniform mat4  MVP;            // Final modelview and projection camera matrix: CSS -> WCS    
uniform mat4  MVP_inv;        // Final inverse modelview and projection camera matrix: CSS -> WCS 
uniform mat4  P_inv;          // Final inverse camera projection matrix: ECS -> CSS
uniform float R_wcs;          // Rmax: maximum sampling distance (in WCS units)
uniform float factor;         // GI contribution multiplier 
uniform vec3 bbox_min;        // Bounding box limits of the radiance hints volume
uniform vec3 bbox_max;        // 
uniform sampler3D points[4];  // Radiance Hint buffers:
                              // points[0]: dist_min, dist_max, dist_ave, 1.0
                              // points[1]: Lr(1,1) Lr(1,-1) Lr(1,0) Lr(0,0) 
                              // points[2]: Lg(1,1) Lg(1,-1) Lg(1,0) Lg(0,0) 
                              // points[3]: Lb(1,1) Lb(1,-1) Lb(1,0) Lb(0,0) 

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

vec3 VectorECS2WCS(in vec3 sample) 
{ 
    vec4 vector_WCS = transpose(P_inv*MVP)*vec4(sample,0);
    return vector_WCS.xyz; 
} 

vec3 PointCSS2WCS(in vec3 sample) 
{ 
    vec4 p_wcs = MVP_inv*vec4(sample,1.0); 
    return p_wcs.xyz/p_wcs.w; 
} 

// ----------------  Main shader function -----------------------------------------

void main(void) 
{ 
    // Accumulated global illumination, initialized at (0,0,0)
    vec3 GI = vec3(0.0,0.0,0.0); 
    
    // Blending factor for current GI estimation frame
    // If blending < 1, current result is blended with previous GI values (temporal smoothing) 
    float blending = 0.8; 

    // Discard background fragments
    float depth = texture2D(RT_depth,gl_TexCoord[0].st).r;
    if (depth==1.0) 
    {
	gl_FragColor=vec4(0.0,0.0,0.0,1.0);//discard; 
	return; 
    } 
	
    vec3 extents = bbox_max-bbox_min; 
    ivec3 sz = textureSize(points[0],0); 

    // convert fragment position and normal to WCS
    vec3 pos_css = vec3(2.0*gl_TexCoord[0].x-1.0, 2.0*gl_TexCoord[0].y-1.0, 2*depth-1.0);
    vec3 pos_wcs = PointCSS2WCS(pos_css); 
    vec4 ntex = texture2D(RT_normals,gl_TexCoord[0].st); 
    vec3 normal_ecs = ntex.xyz; 
    normal_ecs.x = normal_ecs.x*2.0-1.0; 
    normal_ecs.y = normal_ecs.y*2.0-1.0; 
    normal_ecs = normalize(normal_ecs); 
    vec3 normal_wcs = normalize(VectorECS2WCS(normal_ecs)); 
	
    // determine volume texture coordinate of current fragment location
    vec3 uvw = (pos_wcs-bbox_min)/extents; 
    
    float denom = 0.05; 

    // Sample the RH volume at 4 locations, one directly above the shaded point,
    // three on a ring 80degs away from the normal direction. All samples are
    // moved away from the shaded point to avoid sampling RHs behind the surface. 
    //    
    // You can introduce a random rotation of the samples around the normal:
    vec3 rnd = vec3(0,0,0); //texture3D(Noise,25*uvw).xyz - vec3(0.5,0.5,0.5); 

    // Generate the sample locations;
    vec3 v_rand = vec3(0.5, 0.5, 0.5); 
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
        
        vec3 rh_pos    = bbox_min+extents*uvw_new; 
	vec3 path = rh_pos - pos_wcs; 
	float dist = length(path); 
	float rel_dist = dist/R_wcs; 
	dist/=length(extents/sz); 
	path = normalize(path); 
	float contrib = dist>0.005?1.0:0.0; 
	GI+= contrib*SH2RGB (rh_shr, rh_shg, rh_shb, -normal_wcs); 
	denom+=contrib; 
    } 
    GI*=factor/denom; 
    
    // Note: tone mapping is normally required on the final GI buffer
    gl_FragColor = vec4(GI,blending); 
}
