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
using Eigen::Vector4f;

using namespace std;

// Constants based on what was given in the task
static const int keyframe_count = 5;
static const string in_directory = "keyframes/";
static const string out_directory = "output/";
static const string frame_name = "bunny";
static const int keyframe_nums[keyframe_count] = {0, 5, 10, 15, 20};

struct Point 
{
    float x;
    float y;
    float z;
};

struct Face
{
    int p1;
    int p2;
    int p3;
};

struct Frame
{
    int number;
    vector<Point *> points;
};


// Vectors used to store keyframes, interpolated frames, and the faces 
vector<Frame *> keyframes;
vector<Frame *> interpolated_frames;
vector<Face *> faces;


// Prints usage message if user calls the program with incorrect arguments
void usage(string filename) {
    cerr << "usage: " << filename << endl;
    exit(1);
}


// Takes any integer [0, 99] and returns it as a 2-digit string
string intTo2DigitString(int num) {
    if (num < 10) {
        return string("0" + to_string(num));
    }
    return to_string(num);
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

void test_help(int idx) {
    // Makes keyframe
    Frame *keyframe = new Frame;
    keyframe->number = keyframe_nums[idx];

    // Opens file
    /* WHEN I REMOVE OPENING AND CLOSING THE FILE, THERE'S NO ERROR, 
     * BUT ONCE I ADD IN THE FILE, IT ERRORS ON THE SECOND ITERATION, 
     * WHICH DOESNT MAKE SENSE TO ME BECAUSE I DONT EVEN READ ANYTHING 
     * FROM THE FILE
     */
    string path = in_directory + frame_name + intTo2DigitString(keyframe_nums[idx]) + ".obj";
    ifstream obj_file(path.c_str());
    if (obj_file.fail()) {
        throw invalid_argument("Could not read input obj file '" + path + "'.");
    }
 
    // adds a singular point that isn't even read in from the file
    {
        Point *p = new Point;
        p->x = 10;
        p->y = 20;
        p->z = 30;
        cerr << idx << " keyframe addr: " << keyframe 
                << " point vector addr: " << &(keyframe->points) 
                << " point addr: " << p << endl;

        /* STILL ERRORS RIGHT HERE ON THE SECOND ITERATION.
         * FOR SOME REASON THE FIRST ITERATION READS IN THE POINT FINE 
         * BOTH OPENING AND CLOSING THE FILE, BUT ONCE IT GETS TO THE SECOND
         * ITERATION IT ERRORS ON THIS LINE
         */
        keyframe->points.push_back(p);
    }

    // Doesnt even add the keyframe to our vector 

    // Closes file
    obj_file.close();
}
void test() {
    for (int idx = 0; idx < keyframe_count; idx++) {
        test_help(idx);
    }
}



void read_in_keyframe(int kf_idx, string obj_file_path) {
    Frame *keyframe = new Frame;
    keyframe->number = keyframe_nums[kf_idx];

    // Opens file for reading
    string buffer;
    ifstream obj_file(obj_file_path.c_str());
    if (obj_file.fail()) {
        throw invalid_argument("Could not read input obj file '" + obj_file_path + "'.");
    }

    // Stores intermediate input while reading in parameters
    vector<string> line;

    // Reads in each frame
    while (getline(obj_file, buffer)) {
        line.clear();
        splitBySpace(buffer, line);

        // Reads in vertices
        if (line[0][0] == 'v') {
            Point *p = new Point;
            p->x = stof(line[1]);
            p->y = stof(line[2]);
            p->z = stof(line[3]);
            keyframe->points.push_back(p);
        }

        // Reads in faces but only on the first loop iteration
        if (line[0][0] == 'f') {
            if (kf_idx != 0) {
                break;
            }

            Face *f = new Face;
            f->p1 = stoi(line[1]);
            f->p2 = stoi(line[2]);
            f->p3 = stoi(line[3]);
            faces.push_back(f);
        }
    }

    // Adds the keyframe to our vector and closes the file
    keyframes.push_back(keyframe);
    obj_file.close();
}


// Reads in the keyframes 
void read_input(void) 
{
    for (int idx = 0; idx < keyframe_count; idx++) {
        string path = in_directory + frame_name + intTo2DigitString(keyframe_nums[idx]) + ".obj";
        read_in_keyframe(idx, path);
    }
}


// Uses Catmull-Rom Spline to interpolate a component given values for p-1, p, p+1, & p+2
float interpolate_component(Vector4f vec_u, float comp_m1, float comp, float comp_a1, float comp_a2) 
{
    // Basis Matrix used with splines to interpolate components 
    Eigen::Matrix4f mat_B;
    mat_B << 0, 2, 0, 0,
        -1, 0, 1, 0,
        2, -5, 4, -1,
        -1, 3, -3, 1;
    mat_B *= 0.5;

    // Calculates vector p given the componenets for p-1, p, p+1, & p+2
    Vector4f vec_p;
    vec_p << comp_m1, comp, comp_a1, comp_a2;
    
    // Returns the resulting interpolated component f(u) = u dot B * p
    return vec_u.dot(mat_B * vec_p);
}


// Interpolates the frames between 2 keyframes indicated by idx_p and idx_pa1
void interpolate_gap(int idx_pm1, int idx_p, int idx_pa1, int idx_pa2) {
    // Gets the keyframes for p-1, p, p+1, & p+2
    Frame *fm1 = keyframes[idx_pm1];
    Frame *f = keyframes[idx_p];
    Frame *fa1 = keyframes[idx_pa1];
    Frame *fa2 = keyframes[idx_pa2];

    // Interpolates a frame for every # between the idx_p and idx_pa1 keyframes
    for (int frame_num = f->number + 1; frame_num < fa1->number; frame_num++) {
        Frame* f = new Frame;
        f->number = frame_num;
        
        // Calculates vector u as the input to the Catmull-Rom Spline function f(u)
        float u = (frame_num - f->number) * 1.0 / (fa1->number - f->number);
        Vector4f vec_u;
        vec_u << 1, u, u*u, u*u*u;

        // Interpolates every point for the frame
        for (int point_idx = 0; point_idx < f->points.size(); point_idx++) {
            Point *p = new Point;

            p->x = interpolate_component(vec_u, fm1->points[point_idx]->x,
                                              f->points[point_idx]->x,
                                              fa1->points[point_idx]->x,
                                              fa2->points[point_idx]->x);

            p->y = interpolate_component(vec_u, fm1->points[point_idx]->y,
                                              f->points[point_idx]->y,
                                              fa1->points[point_idx]->y,
                                              fa2->points[point_idx]->y);
            
            p->z = interpolate_component(vec_u, fm1->points[point_idx]->z,
                                              f->points[point_idx]->z,
                                              fa1->points[point_idx]->z,
                                              fa2->points[point_idx]->z);

            f->points.push_back(p);
        }
    }
}


// Gets the previous frame index with a min of 0
int prev_keyframe_idx(int idx) 
{
    return (idx == 0) ? 0 : idx - 1;
}


// Gets the next frame index with a max of keyframes.size() - 1
int next_keyframe_idx(int idx) 
{
    return (idx + 1 == keyframes.size()) ? idx : idx + 1;
}


// Interpolates all the frames
void interpolate_frames(void)
{
    for (int idx = 0; idx < keyframe_count; idx++) {
        int idx_pm1 = prev_keyframe_idx(idx);
        int idx_pa1 = next_keyframe_idx(idx);
        int idx_pa2 = next_keyframe_idx(idx_pa1);
        interpolate_gap(idx_pm1, idx, idx_pa1, idx_pa2);
    }
}


// Output all the interpolated frames
void output_interpolated_frames(void) 
{
    for (int frame_idx = 0; frame_idx < interpolated_frames.size(); frame_idx++) {
        // Gets the interpolated frame for every interpolated frame
        Frame *frame = interpolated_frames[frame_idx];
        string obj_file_path = out_directory + frame_name + intTo2DigitString(frame->number) + ".obj";

        // Creates a file in the out directory for the frame
        string buffer;
        ofstream obj_file(obj_file_path.c_str());
        if (!obj_file.is_open()) {
            throw invalid_argument("Could not open output obj file '" + obj_file_path + "'.");
        } 

        // Writes all the interpolated points
        for (int point_idx = 0; point_idx < frame->points.size(); point_idx++) {
            Point *p = frame->points[point_idx];
            obj_file << "v " << p->x << " " << p->y << " " << p->z << endl;
        }
        
        // Writes all the faces
        for (int face_idx = 0; face_idx < faces.size(); face_idx++) {
            Face *f = faces[face_idx];
            obj_file << "f " << f->p1 << " " << f->p2 << " " << f->p3 << endl;
        }

        obj_file.close();
    }
}


// Frees all heap-allocated data for the program
void destruct() {
    // Frees each keyframe along with all of its points
    for (int frame_idx = 0; frame_idx < keyframe_count; frame_idx++) {
        Frame *frame = keyframes[frame_idx];
        for (int point_idx = 0; point_idx < frame->points.size(); point_idx++) {
            delete frame->points[point_idx];
        }
        delete frame;
    }


    // Frees each interpolated frame along with all of its points
    for (int frame_idx = 0; frame_idx < interpolated_frames.size(); frame_idx++) {
        Frame *frame = interpolated_frames[frame_idx];
        for (int point_idx = 0; point_idx < frame->points.size(); point_idx++) {
            delete frame->points[point_idx];
        }
        delete frame;
    }


    // Frees every face
    for (int face_idx = 0; face_idx < faces.size(); face_idx++) {
        delete faces[face_idx];
    }
}


/* Comments for the individual functions of main are above */
int main(int argc, char* argv[])
{
    // Code added for debugging
    bool debug_testing = true;

    if (debug_testing) {
        test();
    } else {

    // Main program
    if (argc != 1) {
        usage(argv[0]);
    }

    read_input();

    interpolate_frames();

    output_interpolated_frames();

    destruct();

    }
}
