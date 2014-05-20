uniform sampler3D SHR;
uniform sampler3D SHG;
uniform sampler3D SHB;

in vec3 uvw;
out vec4 FragColor;

void main()
{
    float b = texture(SHR, uvw).w;
    float g = texture(SHG, uvw).w;
    float r = texture(SHB, uvw).w;
    FragColor = vec4(r, g, b, 1.0);
}
