#ifndef GL_ARB_bindless_texture
uniform sampler2D tex;
#endif

#ifdef GL_ARB_bindless_texture
flat in uvec2 handle;
#endif
in vec2 uv;
in vec4 color;
out vec4 FragColor;

void main(void)
{
#ifdef GL_ARB_bindless_texture
    vec4 col = texture(sampler2D(handle), uv);
#else
    vec4 col = texture(tex, uv);
#endif
    if (col.a < 0.5) discard;
    FragColor = vec4(1.);
}
