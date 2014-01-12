#pragma once

#include "ofMain.h"
#include "ofxTriangleMesh.h"
#include "ofxRayTriangleIntersection.h"

enum STATE{
    TO_CREATATION   = 0,
    CREATATAING     = 1,
    TO_PAINT        = 2,
    PAINTING        = 3,
    TO_CUT          = 4,
    CUTTING         = 5,
    TO_DRAW_RING    = 6,
    DRAWING_RING    = 7,
    TO_REMOVE_TRI   = 8,
    TO_EXTRUDE      = 9,
    
    ROTATE          = 15,
    TRANSLATE       = 16
};

enum COLOR{
    NONE,
    WHITE,
    GRAY,
    BLACK,
    YELLOW,
};

class of_edge{

public:
    of_edge(){}
    of_edge(ofPoint pt1, ofPoint pt2) { p[0] = pt1; p[1] = pt2; }
    of_edge(ofPoint input[2]) { p[0]= input[0];p[1]= input[1]; }
    
    ~of_edge() { p[0].~ofVec3f();p[1].~ofVec3f(); }
    
    // public function
    void draw() { ofLine(p[0], p[1]); }
    
    // class member
    ofPoint p[2];
};

class of_triangle{

public:
    of_triangle() {}
    of_triangle(ofPoint input[3]) { p[0]= input[0];p[1]= input[1];p[2]= input[2];chordal_axis.clear(); }
    of_triangle(of_triangle const &input)
    {
        p[0]= input.p[0];p[1]= input.p[1];p[2]= input.p[2];
        type = input.type;
        counter[0] = input.counter[0];counter[1] = input.counter[1];counter[2] = input.counter[2];
        chordal_axis = input.chordal_axis;
        line_seg = input.line_seg;
        normal[0]=input.normal[0];normal[1]=input.normal[1];normal[2]=input.normal[2];
    }
    of_triangle(ofPoint p0, ofPoint p1, ofPoint p2){
        p[0] = p0;
        p[1] = p1;
        p[2] = p2;
    }

    ~of_triangle(){ p[0].~ofVec3f();p[1].~ofVec3f();p[2].~ofVec3f(); }
    
    // public function
    void draw_triangle() { ofTriangle(p[0], p[1], p[2]); }
    void draw_wireframe(){ ofLine(p[0], p[1]);ofLine(p[1], p[2]);ofLine(p[2], p[0]); }
    void copyNormal(ofVec3f norm[3]){
        normal[0] = norm[0];
        normal[1] = norm[1];
        normal[2] = norm[2];
    }
    
    // class member
    ofPoint p[3];
    int type;
    int counter[3]={0,0,0};
    vector<of_edge> chordal_axis,line_seg;
    ofVec3f normal[3];
    
    int numSplit = -1;
    COLOR vColor[3] = {NONE, NONE, NONE};
    COLOR tColor = NONE;
};

class testApp : public ofBaseApp{
    
public:
    // About OF
    void setup();
    void update();
    void draw();
    
    void keyPressed  (int key);
    void keyReleased(int key);
    void mouseMoved(int x, int y );
    void mouseDragged(int x, int y, int button);
    void mousePressed(int x, int y, int button);
    void mouseReleased(int x, int y, int button);
    void windowResized(int w, int h);
    void dragEvent(ofDragInfo dragInfo);
    void gotMessage(ofMessage msg);
    
    void prune_1();
    void prune_2();
    bool check_advance(vector<ofPoint>& collect_point,ofPoint& center,float r);
    bool check_same_egde(ofPoint& e1p1,ofPoint& e1p2,ofPoint& e2p1,ofPoint& e2p2);
    void add_fan_triangle_1(vector<of_edge>& collect_edge,ofPoint center);
    void add_fan_triangle_2(ofPoint P1,ofPoint P2,ofPoint P3);
    void add_fan_triangle_3(ofPoint P1,ofPoint P2,ofPoint P3);
    void elevate();
    bool check_notinside(vector<ofPoint>& added_point,ofPoint check_point);
    void cast();
    void clone();
    void quarter_oval();
    void paint_line();
    void cut_plane();
    void cut();
    void cut_construct(ofVec3f& plane_normal);
    
    // operating function
    void rotate(float theta, ofVec3f dir);
    void translate(float dist, ofVec3f dir);
    
    // To createRing
    void sortLineSeg();
    void seperateTri();
    void drawRingTriangle(of_triangle* tri, int baseIdx);
    void simplifyLine();
    vector<of_edge> simplifyOneLine(of_triangle* tri, int first, int last);
    void createRing();
    
    // class member
    ofPolyline line;
    
    ofxTriangleMesh mesh;
    
    vector<of_triangle> Tlist;
    int T_num=0;
    
    ofLight pointLight;
    ofCamera cam;
    ofMaterial material;
    bool release = false;
    bool draw_c = false;
    bool elevated_T = false;
    bool enableFace = true;
    
    STATE state = TO_CREATATION;
    STATE stateSave = TO_CREATATION;
    
    int old_x = 0;
    int old_y = 0;
    
    // for extrusion
    vector<int> triBelongToRing;
    vector<int> triInsideRing;
    
    COLOR colorInsideRing = NONE;
    
    // for test
    bool isTest = false;
};
