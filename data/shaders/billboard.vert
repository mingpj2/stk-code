uniform mat4 ModelViewMatrix;
uniform vec3 Position;
uniform vec2 Size;

#if __VERSION__ >= 330
layout(location = 0) in vec2 Corner;
layout(location = 3) in vec2 Texcoord;
#else
in vec2 Corner;
in vec2 Texcoord;
#endif

out vec2 uv;

void main(void)
{
    uv = Texcoord;
    vec4 Center = ModelViewMatrix * vec4(Position, 1.);
    gl_Position = ProjectionMatrix * (Center + vec4(Size * Corner, 0., 0.));
}
