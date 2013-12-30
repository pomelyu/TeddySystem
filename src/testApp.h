#pragma once

#include "ofMain.h"
#include "ofxTriangleMesh.h"

class of_edge
{
public:
    of_edge(){}
    of_edge(ofPoint input[2]) {p[0]= input[0];p[1]= input[1];}
    ~of_edge(){p[0].~ofVec3f();p[1].~ofVec3f();}
    void draw() {ofLine(p[0], p[1]);}
    ofPoint p[2];
};

class of_triangle
{
public:
    of_triangle(ofPoint input[3]) {p[0]= input[0];p[1]= input[1];p[2]= input[2];chordal_axis.clear();}
    ~of_triangle(){p[0].~ofVec3f();p[1].~ofVec3f();p[2].~ofVec3f();}
    void draw_triangle() {ofTriangle(p[0], p[1], p[2]);}
    void draw_wireframe(){ofLine(p[0], p[1]);ofLine(p[1], p[2]);ofLine(p[2], p[0]);}
    ofPoint p[3];
    int type;
    int counter[3]={0,0,0};
    vector<of_edge> chordal_axis;
    ofVec3f normal[3];
};

class testApp : public ofBaseApp{
    
public:
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
    
    
    
};
