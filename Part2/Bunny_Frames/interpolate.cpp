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
static const string keyframe_nums[keyframe_count] = {"00", "05", "10", "15", "20"};

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
    vector<Point> points;
};


// Vectors used to store keyframes, interpolated frames, and the faces 
vector<Frame *> keyframes;
vector<Frame *> interpolated_frames;
vector<Face> faces;


// Prints usage message if user calls the program with incorrect arguments
void usage(string filename) {
    cerr << "usage: " << filename << endl;
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


// Reads in the keyframes 
void read_keyframes(void) 
{
    // Mallocs a new frame for each key frame
    for (int idx = 0; idx < keyframe_count; idx++) {
        Frame *keyframe = (Frame*)malloc(sizeof(Frame));
        keyframe->number = stoi(keyframe_nums[idx]);

        // Opens the input file for each keyframe
        string obj_file_path = in_directory + frame_name + keyframe_nums[idx] + ".obj";
        string buffer;
        ifstream obj_file;
        obj_file.open(obj_file_path.c_str(), ifstream::in);
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
                cerr << buffer << endl;
                Point p;
                p.x = stof(line[1]);
                p.y = stof(line[2]);
                p.z = stof(line[3]);
                keyframe->points.push_back(p);
            }

            // Reads in faces but only on the first loop iteration
            if (line[0][0] == 'f') {
                if (idx != 0) {
                    cerr << buffer << "Breaking.." << endl;

                    break;
                }

                cerr << buffer << endl;

                Face f;
                f.p1 = stoi(line[1]);
                f.p2 = stoi(line[2]);
                f.p3 = stoi(line[3]);
                faces.push_back(f);
            }
        }

        // Adds the keyframe to our vector and closes the file input once the frame is read in
        keyframes.push_back(keyframe);
        obj_file.close();
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
        Frame* f = (Frame *)malloc(sizeof(Frame));
        f->number = frame_num;

        // Calculates vector u as the input to the Catmull-Rom Spline function f(u)
        float u = (frame_num - f->number) * 1.0 / (fa1->number - f->number);
        Vector4f vec_u;
        vec_u << 1, u, u*u, u*u*u;

        // Interpolates every point for the frame
        for (int point_idx = 0; point_idx < f->points.size(); point_idx++) {
            Point p;

            p.x = interpolate_component(vec_u, fm1->points[point_idx].x,
                                              f->points[point_idx].x,
                                              fa1->points[point_idx].x,
                                              fa2->points[point_idx].x);

            p.y = interpolate_component(vec_u, fm1->points[point_idx].y,
                                              f->points[point_idx].y,
                                              fa1->points[point_idx].y,
                                              fa2->points[point_idx].y);
            
            p.z = interpolate_component(vec_u, fm1->points[point_idx].z,
                                              f->points[point_idx].z,
                                              fa1->points[point_idx].z,
                                              fa2->points[point_idx].z);

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


// Takes any integer [0, 99] and returns it as a 2-digit string
string intTo2DigitString(int num) {
    if (num < 10) {
        return string("0" + to_string(num));
    }
    return to_string(num);
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
            Point p = frame->points[point_idx];
            obj_file << "v " << p.x << " " << p.y << " " << p.z << endl;
        }
        
        // Writes all the faces
        for (int face_idx = 0; face_idx < faces.size(); face_idx++) {
            Face f = faces[face_idx];
            obj_file << "f " << f.p1 << " " << f.p2 << " " << f.p3 << endl;
        }

        obj_file.close();
    }
}


// Frees all heap-allocated data for the program
void destruct() {
    for (int idx = 0; idx < keyframe_count; idx++) {
        free(keyframes[idx]);
    }

    for (int idx = 0; idx < interpolated_frames.size(); idx++) {
        free(interpolated_frames[idx]);
    }
}


/* Comments for the individual functions of main are above */
int main(int argc, char* argv[])
{
    if (argc != 1) {
        usage(argv[0]);
    }

    read_keyframes();

    interpolate_frames();

    output_interpolated_frames();

    destruct();
}
