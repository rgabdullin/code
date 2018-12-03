#include <mgl2/glut.h>

int main()
{
    mglGLUT gr;
    gr.FPlot("sin(pi*x)");
    gr.WriteFrame("test.png");
    return 0;
}