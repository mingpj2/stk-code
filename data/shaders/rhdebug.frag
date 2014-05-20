uniform sampler3D SHR;
uniform sampler3D SHG;
uniform sampler3D SHB;

in vec3 uvw;
out vec4 FragColor;

void main()
{
    float r = texture(SHR, uvw).x;
    float g = texture(SHG, uvw).x;
    float b = texture(SHB, uvw).x;
    FragColor = vec4(r, g, b, 1.0);
}
