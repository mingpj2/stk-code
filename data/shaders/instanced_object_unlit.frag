#ifndef GL_ARB_bindless_texture
uniform sampler2D tex;
#endif

#ifdef GL_ARB_bindless_texture
flat in sampler2D handle;
#endif

in vec2 uv;
in vec4 color;
out vec4 FragColor;

void main(void)
{
#ifdef GL_ARB_bindless_texture
    vec4 col = texture(handle, uv);
#ifdef SRGBBindlessFix
    col.xyz = pow(col.xyz, vec3(2.2));
#endif
#else
    vec4 col = texture(tex, uv);
#endif
    col.xyz *= pow(color.xyz, vec3(2.2));
    if (col.a < 0.5) discard;
    FragColor = vec4(col.xyz, 1.);
}
