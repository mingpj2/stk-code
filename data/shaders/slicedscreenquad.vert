in vec2 Position;
in vec2 Texcoord;
uniform int slice;

#if __VERSION__ >= 130
out vec2 uv;
#else
varying vec2 uv;
#endif


void main() {
    gl_Layer = slice;
	uv = Texcoord;
	gl_Position = vec4(Position, 0., 1.);
}
