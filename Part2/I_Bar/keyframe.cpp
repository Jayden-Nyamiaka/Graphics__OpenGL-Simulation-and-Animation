/* Libraries neeeded for OpenGL*/
#include <GL/glut.h>

/* Libraires neede for math computations */
#include <math.h>
#define _USE_MATH_DEFINES

/* Standard libraries that are just generally useful. */
#include <iostream>
#include <vector>

/* Libraries used for file parsing */
#include <fstream>
#include <sstream>

/* Eigen Library included for ArcBall implementation */
#include <Eigen/Dense>
using Eigen::Vector3f;
using Eigen::Matrix4f;

using namespace std;

struct Quarternion
{
    float real;
    float im[3];
};

struct Frame
{
    int number;
    float translation[3];
    float scaling[3];
    Quarternion rotation;
};

/*Quarternion getIdentityQuarternion(void) {
    Quarternion q;
    q.real = 1.0f;
    q.im[0] = 0.0f;
    q.im[1] = 0.0f;
    q.im[2] = 0.0f;
    return q;
}*/

/* This is a code snippet for drawing the "I-bar". The 'quadratic' object can be made
 * into a global variable in your program if you want. Line 13, where 'quadratic' gets
 * initialized should be done in your 'init' function or somewhere close to the start
 * of your program.
 *
 * The I-bar is an object that Prof. Al Barr once used in one of his papers to demonstrate
 * an interpolation technique. You might call it "Al Barr's I-Bar" ;)
 */

/* Needed to draw the cylinders using glu */
GLUquadricObj *quadratic;

///////////////////////////////////////////////////////////////////////////////////////////////////

/* Hardcoded camera and frustum parameters */
float cam_position[3] = {0.0f, 0.0f, 40.0f};
// No camera rotation - no orientation axis or angle
float near_param = 1.0f,
    far_param = 60.0f,
    left_param = -1.0f,
    right_param = 1.0f,
    bottom_param = -1.0f,
    top_param = 1.0f;

///////////////////////////////////////////////////////////////////////////////////////////////////

/* Frame parameters read in from the input .script file */
int frame_num;
int frame_count;
vector<Frame*> frames;

///////////////////////////////////////////////////////////////////////////////////////////////////

/* Parameters for drawing the cylinders */
float cyRad = 0.2;
float cyHeight = 1.0;
int quadStacks = 4;
int quadSlices = 4;

///////////////////////////////////////////////////////////////////////////////////////////////////


void usage(string filename) {
    cerr << "usage: " << filename << " input_file.script xres yres\n\t"
            "xres, yres (screen resolution) must be positive integers\n";
    exit(1);
}


// Tokenizes a string by spaces returning the tokenized elements in a vector
void splitBySpace(string s, vector<string> &split)
{
    stringstream stream(s);

    string buffer;
    while(getline(stream, buffer, ' ')) {
        if (!buffer.empty()) {
            split.push_back(buffer);
        }
    }
}


// Converts a given angle in degrees to radians
float deg2rad(float angle)
{
    return angle * M_PI / 180.0;
}


void parseScriptFile(string scriptfile)
{
    // Ensures the input script is the 
    if (scriptfile.find(".script") == -1) {
        throw invalid_argument("File " + scriptfile + " needs to be a .script file.");
    }

    string buffer;
    ifstream script;
    script.open(scriptfile.c_str(), ifstream::in);
    if (script.fail()) {
        throw invalid_argument("Could not read input script file '" + scriptfile + "'.");
    }

    // These variables store intermediate input while reading in parameters
    vector<string> line;
    Frame *frame;

    /* Reads in all frame parameters  */
    while (getline(script, buffer)) {
        line.clear();
        splitBySpace(buffer, line);

        if (line.size() == 0) {
            break;
        } else if (line.size() == 1) {
            frame_count = stoi(line[0]);
        } else if (line.size() == 2) {
            frame = (Frame *)malloc(sizeof(Frame));
            frame->number = stoi(line[1]);
        } else if (line[0][0] == 't') {
            frame->translation[0] = stof(line[1]);
            frame->translation[1] = stof(line[2]);
            frame->translation[2] = stof(line[3]);
        } else if (line[0][0] == 's') {
            frame->scaling[0] = stof(line[1]);
            frame->scaling[1] = stof(line[2]);
            frame->scaling[2] = stof(line[3]);
        } else if (line[0][0] == 'r') {
            float rot_x = stof(line[1]);
            float rot_y = stof(line[2]);
            float rot_z = stof(line[3]);
            float theta = deg2rad( stof(line[4]) );
            double sinHalfTheta = sin(0.5 * theta);

            frame->rotation.real = cos(0.5 * theta);
            frame->rotation.im[0] = rot_x * sinHalfTheta;
            frame->rotation.im[1] = rot_y * sinHalfTheta;
            frame->rotation.im[2] = rot_z * sinHalfTheta;

            frames.push_back(frame);
        }
    }
    script.close();
}


void init(int &argc, char* argv[])
{
    /* Checks that the user inputted the right parameters into the command line
     * and stores xres, yres, and filename to their respective fields
     */
    string filename = argv[0];
    if (argc != 4) {
        usage(filename);
    }
    string script_file = argv[1];
    int xres = stoi(argv[2]);
    int yres = stoi(argv[3]);
    if (xres <= 0 || yres <= 0) {
        usage(filename);
    }

    // Intializes the GLUT library
    glutInit(&argc, argv);

    // Tells OpenGL we need a double buffer, a RGB pixel buffer, and a depth buffer
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);

    // Creates the OpenGL program window
    glutInitWindowSize(xres, yres);

    // Set the origin of the program window to the top-left corner i.e (0, 0)
    glutInitWindowPosition(0, 0);

    // Names the program window and creates it
    glutCreateWindow("I-bar Animation");

    // Extracts all information from script file entered in command line
    parseScriptFile(script_file);

    // Declares quadratic needed to draw I-bar later
    quadratic = gluNewQuadric();

    // Initializes Rotation Quarternions
    //last_rotation = getIdentityQuarternion();
    //curr_rotation = getIdentityQuarternion();

    // Specifies Gouraud Shading as the shading mode
    glShadeModel(GL_SMOOTH);
    
    // Uses Backface Culling as an optimization when rendering
    glEnable(GL_CULL_FACE);
    glCullFace(GL_BACK);
    
    // Uses Depth Buffering as an optimization when rendering
    glEnable(GL_DEPTH_TEST);
    
     /* Tells OpenGL to automatically normalize our normal
     * vectors before it passes them into the normal arrays discussed below
     */
    glEnable(GL_NORMALIZE);
    
    /* The following two lines tell OpenGL to enable its "vertex array" and
     * "normal array" functionality.
     */
    //glEnableClientState(GL_VERTEX_ARRAY);
    //glEnableClientState(GL_NORMAL_ARRAY);
    
    /* The next 4 lines work with OpenGL's two main matrices: the "Projection
     * Matrix" and the "Modelview Matrix". Only one of these two main matrices
     * can be modified at any given time. We specify the main matrix that we
     * want to modify with the 'glMatrixMode' function.
     *
     * The Projection Matrix is the matrix that OpenGL applies to points in
     * camera space. For our purposes, we want the Projection Matrix to be
     * the perspective projection matrix, since we want to convert points into
     * NDC after they are in camera space.
     *
     * The line of code below:
     */
    glMatrixMode(GL_PROJECTION);
    /* ^tells OpenGL that we are going to modify the Projection Matrix. From
     * this point on, any matrix comamnds we give OpenGL will affect the
     * Projection Matrix. For instance, the line of code below:
     */
    glLoadIdentity();
    /* ^tells OpenGL to set the current main matrix (which is the Projection
     * Matrix right now) to the identity matrix. Then, the next line of code:
     */
    glFrustum(left_param, right_param,
              bottom_param, top_param,
              near_param, far_param);
    /* ^ tells OpenGL to create a perspective projection matrix using the
     * given frustum parameters. OpenGL then post-multiplies the current main
     * matrix (the Projection Matrix) with the created matrix. i.e. let 'P'
     * be our Projection Matrix and 'F' be the matrix created by 'glFrustum'.
     * Then, after 'F' is created, OpenGL performs the following operation:
     *
     * P = P * F
     * 
     * Since we had set the Projection Matrix to the identity matrix before the
     * call to 'glFrustum', the above multiplication results in the Projection
     * Matrix being the perspective projection matrix, which is what we want.
     */
    
    /* The Modelview Matrix is the matrix that OpenGL applies to untransformed
     * points in world space. OpenGL applies the Modelview Matrix to points
     * BEFORE it applies the Projection Matrix.
     * 
     * Thus, for our purposes, we want the Modelview Matrix to be the overall
     * transformation matrix that we apply to points in world space before
     * applying the perspective projection matrix. This means we would need to
     * factor in all the individual object transformations and the camera
     * transformations into the Modelview Matrix.
     *
     * The following line of code tells OpenGL that we are going to modify the
     * Modelview Matrix. From this point on, any matrix commands we give OpenGL
     * will affect the Modelview Matrix.
     *
     * We generally modify the Modelview Matrix in the 'display' function,
     * right before we tell OpenGL to render anything. See the 'display'
     * for details.
     */
    glMatrixMode(GL_MODELVIEW);
    
    /* The next line calls our function that tells OpenGL to initialize some
     * lights to represent our Point Light structs. Further details will be
     * given in the function itself.
     *
     * The reason we have this procedure as a separate function is to make
     * the code more organized.
     */
    //init_lights();
}


void drawIBar()
{   
    glPushMatrix();
    glColor3f(0, 0, 1);
    glTranslatef(0, cyHeight, 0);
    glRotatef(90, 1, 0, 0);
    gluCylinder(quadratic, cyRad, cyRad, 2.0 * cyHeight, quadSlices, quadStacks);
    glPopMatrix();
    
    glPushMatrix();
    glColor3f(0, 1, 1);
    glTranslatef(0, cyHeight, 0);
    glRotatef(90, 0, 1, 0);
    gluCylinder(quadratic, cyRad, cyRad, cyHeight, quadSlices, quadStacks);
    glPopMatrix();
    
    glPushMatrix();
    glColor3f(1, 0, 1);
    glTranslatef(0, cyHeight, 0);
    glRotatef(-90, 0, 1, 0);
    gluCylinder(quadratic, cyRad, cyRad, cyHeight, quadSlices, quadStacks);
    glPopMatrix();
    
    glPushMatrix();
    glColor3f(1, 1, 0);
    glTranslatef(0, -cyHeight, 0);
    glRotatef(-90, 0, 1, 0);
    gluCylinder(quadratic, cyRad, cyRad, cyHeight, quadSlices, quadStacks);
    glPopMatrix();
    
    glPushMatrix();
    glColor3f(0, 1, 0);
    glTranslatef(0, -cyHeight, 0);
    glRotatef(90, 0, 1, 0);
    gluCylinder(quadratic, cyRad, cyRad, cyHeight, quadSlices, quadStacks);
    glPopMatrix();
}


void display(void)
{
    // Clears our Color and Depth Buffers
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    
    // Resets our ModelView Matrix
    glLoadIdentity();

    // Translates by our camera's position
    glTranslatef(-cam_position[0], -cam_position[1], -cam_position[2]);


    /* Our next step is to set up all the lights in their specified positions.
     * Our helper function, 'set_lights' does this for us. See the function
     * for more details.
     *
     * The reason we have this procedure as a separate function is to make
     * the code more organized.
     */
    //set_lights();

    // Draws the I-bar where and how we specified
    drawIBar();
    
    // Swaps our double buffer
    glutSwapBuffers();
}


void reshape(int width, int height)
{
    /* The following two lines of code prevent the width and height of the
     * window from ever becoming 0 to prevent divide by 0 errors later.
     * Typically, we let 1x1 square pixel be the smallest size for the window.
     */
    height = (height == 0) ? 1 : height;
    width = (width == 0) ? 1 : width;
    
    /* The 'glViewport' function tells OpenGL to determine how to convert from
     * NDC to screen coordinates given the dimensions of the window. The
     * parameters for 'glViewport' are (in the following order):
     *
     * - int x: x-coordinate of the lower-left corner of the window in pixels
     * - int y: y-coordinate of the lower-left corner of the window in pixels
     * - int width: width of the window
     * - int height: height of the window
     *
     * We typically just let the lower-left corner be (0,0).
     *
     * After 'glViewport' is called, OpenGL will automatically know how to
     * convert all our points from NDC to screen coordinates when it tries
     * to render them.
     */
    glViewport(0, 0, width, height);
    
    // Rerenders everything that was on the screen
    glutPostRedisplay();
}


// TODO: Later once things renders
void transform_next_frame(void)
{
    if (frame_num == frame_count - 1) {
        // set frame to 0 frame
        return;
    }

    // else just set next frame
}


// Display the next frame if any key is pressed
void key_pressed(unsigned char key, int x, int y)
{
    transform_next_frame();
}


void destruct() {
    for (int i = 0; i < frames.size(); i++) {
        free(frames[i]);
    }
}


int main(int argc, char* argv[])
{
    // Sets up everything we need for the program
    init(argc, argv);

    // Specifies to OpenGL our display function
    glutDisplayFunc(display);

    // Specifies to OpenGL our reshape function
    glutReshapeFunc(reshape);
    
    // No handling for mouse presses or movement

    // Specifies to OpenGL our function for handling key presses
    glutKeyboardFunc(key_pressed);
    
    // Starts out main loop the program
    glutMainLoop();

    // Frees all heap-allocated data
    destruct();
}
