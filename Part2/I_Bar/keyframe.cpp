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
using Eigen::Vector4f;
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

/* Frame parameters */

// Basis Matrix used with splines to interpolate components 
Matrix4f mat_B;
B << 0, 2, 0, 0,
     -1, 0, 1, 0,
     2, -5, 4, -1,
     -1, 3, -3, 1;
B *= 0.5;

// Parameters read in from the input .script file
int frame_count;
vector<Frame*> frames;

// Parameters that track the current state of the frames and which is being displayed
int frame_num; // the # of the frame that we are currently displaying
int frame_idx; // the idx of the frame in our frames vector that we are displaying or using to interpolate as p_i 
Frame currFrame;

// Index parameters that control interpolation between key frames
int idx_pm1, idx_pa1, idx_pa2; // the frame_idx works as idx_p / p_i

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

    // Sets the current frame to frame 0
    frame_num = 0;
    frame_idx = 0;
    currFrame = *(frames[frame_idx]);
    
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
}


// Applies trnasformation for the current frame
void transformIBar(void)
{
    // Applies the translation first
    glTranslatef(currFrame.translation[0], 
                 currFrame.translation[1], 
                 currFrame.translation[2]);

    // Then applies the rotation
    Quarternion q = currFrame.rotation;
    GLfloat rot[16];

    rot[0] = 1.0f - 2.0f * q.im[1] * q.im[1] - 2.0f * q.im[2] * q.im[2];
    rot[1] = 2.0f * (q.im[0] * q.im[1] - q.im[2] * q.real);
    rot[2] = 2.0f * (q.im[0] * q.im[2] + q.im[1] * q.real);
    rot[3] = 0.0f;

    rot[4] = 2.0f * (q.im[0] * q.im[1] + q.im[2] * q.real);
    rot[5] = 1.0f - 2.0f * q.im[0] * q.im[0] - 2.0f * q.im[2] * q.im[2];
    rot[6] = 2.0f * (q.im[1] * q.im[2] - q.im[0] * q.real);
    rot[7] = 0.0f;

    rot[8] = 2.0f * (q.im[0] * q.im[2] - q.im[1] * q.real);
    rot[9] = 2.0f * (q.im[1] * q.im[2] + q.im[0] * q.real);
    rot[10] = 1.0f - 2.0f * q.im[0] * q.im[0] - 2.0f * q.im[1] * q.im[1];
    rot[11] = 0.0f;
    
    rot[12] = 0.0f;
    rot[13] = 0.0f;
    rot[14] = 0.0f;
    rot[15] = 1.0f;

    glMultMatrixf(rot);

    // Finally applies the scaling
    glScalef(currFrame.scaling[0], 
             currFrame.scaling[1], 
             currFrame.scaling[2]);
}

void drawIBar(void)
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

    // Applies the set I-bar transformation to our ModelView Matrix
    transformIBar();

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


// Uses Catmull-Rom Spline to interpolate a component given values for p-1, p, p+1, & p+2
void interpolateComponent(Vector4f vec_u, float comp_m1, float comp, float comp_a1, float comp_a2) 
{
    // Calculates vector p given the componenets for p-1, p, p+1, & p+2
    Vector4f vec_p;
    vec_p << comp_m1, comp, comp_a1
    
    // returns the resulting interpolated component f(u) = u B p
    return vec_u * mat_B * vec_p;
}


void interpolate(void)
{
    // Gets the frames indicated by frame indices for p-1, p, p+1, & p+2
    Frame *fm1 = frames[idx_pm1];
    Frame *f = frames[frame_idx];
    Frame *fa1 = frames[idx_pa1];
    Frame *fa2 = frames[idx_pa2];

    // Calculates vector u as the input to the Catmull-Rom Spline function f(u)
    float u = (frame_num - f->number) / (fa1->number - f->number);
    Vector4f vec_u;
    vec_u << 1, u, u*u, u*u*u;

    // Sets the frame number
    currFrame.number = frame_num;

    // Interpolates all components

    // Translation
    currFrame.translation[0] = interpolateComponent(vec_u, fm1->translation[0],
                                                           f->translation[0],
                                                           fa1->translation[0],
                                                           fa2->translation[0]);
    currFrame.translation[1] = interpolateComponent(vec_u, fm1->translation[1],
                                                           f->translation[1],
                                                           fa1->translation[1],
                                                           fa2->translation[1]);
    currFrame.translation[2] = interpolateComponent(vec_u, fm1->translation[2],
                                                           f->translation[2],
                                                           fa1->translation[2],
                                                           fa2->translation[2]);

    // Scaling
    currFrame.scaling[0] = interpolateComponent(vec_u, fm1->scaling[0],
                                                       f->scaling[0],
                                                       fa1->scaling[0],
                                                       fa2->scaling[0]);
    currFrame.scaling[1] = interpolateComponent(vec_u, fm1->scaling[1],
                                                       f->scaling[1],
                                                       fa1->scaling[1],
                                                       fa2->scaling[1]);
    currFrame.scaling[2] = interpolateComponent(vec_u, fm1->scaling[2],
                                                       f->scaling[2],
                                                       fa1->scaling[2],
                                                       fa2->scaling[2]);
    
    // Rotation
    Vector3f rotIm;
    rotIm[0] = interpolateComponent(vec_u, fm1->rotation.im[0],
                                           f->rotation.im[0],
                                           fa1->rotation.im[0],
                                           fa2->rotation.im[0]);
    rotIm[1] = interpolateComponent(vec_u, fm1->rotation.im[1],
                                           f->rotation.im[1],
                                           fa1->rotation.im[1],
                                           fa2->rotation.im[1]);
    rotIm[2] = interpolateComponent(vec_u, fm1->rotation.im[2],
                                           f->rotation.im[2],
                                           fa1->rotation.im[2],
                                           fa2->rotation.im[2]);
    rotIm.normalize();
    currFrame.rotation.im[0] = rotIm[0];
    currFrame.rotation.im[1] = rotIm[1];
    currFrame.rotation.im[2] = rotIm[2];

}


int prev_frame_idx(int idx) 
{
    return (idx == 0) ? frames.size() - 1 : idx - 1;
}


int next_frame_idx(int idx) 
{
    return (idx + 1 == frames.size()) ? 0 : idx + 1;
}


void next_frame(void)
{
    // Increments our frame number
    frame_num++;

    // If we've reached the last frame of our animation, go back to the first frame
    if (frame_num == frame_count) {
        frame_num = 0;
        frame_idx = 0;
        currFrame = *(frames[frame_idx]);
    }

    // If the next frame is a keyframe, set our current frame to that keyframe
    else if (frame_num == frames[frame_idx + 1]->number) {
        frame_idx++;
        currFrame = *(frames[frame_idx]);
    }

    // Else the frame is not a keyframe and must be interpolated
    else {
        // if the prev frame displayed is a keyframe, set new interpolation indices
        if (frame_num - 1 == frames[frame_idx]->number) {
            idx_pm1 = prev_frame_idx(frame_idx);
            idx_pa1 = next_frame_idx(frame_idx);
            idx_pa2 = next_frame_idx(idx_pa1);
        }
        interpolate();
    }
}


// Display the next frame if any key is pressed
void key_pressed(unsigned char key, int x, int y)
{
    next_frame();
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