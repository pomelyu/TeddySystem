//
//  extrude.cpp
//  TeddySystem
//
//  Created by Chien Chin-yu on 2014/1/8.
//
//

#include "testApp.h"
#include "testApp.cpp"

#pragma mark -
#pragma mark About extrusion

void testApp::seperateTri(){
    
    // Draw all triangle in ring to WHITE
    for (int i = 0; i < triBelongToRing.size(); i++)
        Tlist[triBelongToRing[i]].tColor =WHITE;
    
    // Draw one of vertex in a ring triangle to BLACK
    // Draw the triangle to GRAY
    // GRAY means the vertex of triangle has been drawed
    of_triangle* tri = &(Tlist[triBelongToRing[0]]);
    float level = tri -> line_seg[0].p[0].y;
    float diff = tri -> p[0].y - level;
    
    tri -> vColor[0] = GRAY;
    for (int i = 1; i < 3; i++) {
        if ((tri -> p[i].y - level) * diff > 0)
            tri -> vColor[i] = GRAY;
    }
    tri -> tColor = GRAY;
    
    vector<of_triangle*> triStack;
    triStack.push_back(tri);
    ofPoint pt = tri -> p[0];
    
    while (triStack.size() > 0) {
        tri = triStack.back();
        int idx = 0;
        for (; idx < 3; idx++) {
            
            // get the first none BLACL vertex
            if (tri -> vColor[idx] == GRAY){
                tri -> vColor[idx] = BLACK;
                pt = tri -> p[idx];
                
                
                for (int i = 0; i < Tlist.size(); i++) {
                    tri = &(Tlist[i]);
                    
                    // triangle in Ring
                    if (tri -> tColor == WHITE) {
                        // check if vertex in the triangle
                        for (int j = 0; j < 3; j++){
                            if (check_point_same(pt, tri -> p[j])) {
                                // draw the three vertex to appropriate color
                                tri -> vColor[j] == BLACK;
                                level = tri -> line_seg[0].p[0].y;
                                diff = tri -> p[j].y - level;
                                for (int k = 0; k < 3; k++) {
                                    if ((tri -> p[k].y - level) * diff > 0)
                                        tri -> vColor[k] = GRAY;
                                }
                                triStack.push_back(tri);
                                break;
                            }
                        }
                    }
                    // adjancent triangle
                    else if (tri -> tColor == NONE){
                        // check if vertex in the triangle
                        for (int j = 0; j < 3; j++) {
                            if (check_point_same(pt, tri -> p[j])) {
                                for (int k = 0; k < 3; k++)
                                    tri -> vColor[k] = GRAY;
                                tri -> vColor[j] = BLACK;
                                tri -> tColor = BLACK;
                                triStack.push_back(tri);
                                break;
                            }
                        }
                    }
                }
                break;
                
            }
        }
        // If the triangle has no GRAY point
        // Pop it
        if (idx == 3)
            triStack.pop_back();
    }
}


int testApp::diffColorTriHasThisPoint(ofPoint pt, COLOR color){
    bool toCheckTri = true;
    int k = -1;
    for (int i = 0; i < Tlist.size(); i++) {
        of_triangle* tri = &(Tlist[i]);
        if ((*tri).tColor == NONE && toCheckTri) {
            for (int j = 0; j < 3; j++) {
                if ((*tri).p[j] == pt) {
                    (*tri).vColor[j] = color;
                    (*tri).tColor = color;
                    k = i;
                    toCheckTri = false;
                }
            }
        }
        else{
            for (int j = 0; j < 3; j++) {
                if ((*tri).p[j] == pt) {
                    (*tri).vColor[j] = color;
                }
            }
        }
    }
    return k;
}


void testApp::createRing(){
    
    // reconstruct the triangle along the ring
    for (int i = 0; i < triBelongToRing.size(); i++) {
        of_triangle tri = Tlist[triBelongToRing[i]];
        float level = tri.line_seg[0].p[0].y;
        float dif[3];
        dif[0] = tri.p[0].y - level;
        dif[1] = tri.p[1].y - level;
        dif[2] = tri.p[2].y - level;
        
        if (dif[0] * dif[1] * dif[2] < 0) { // two vertex in ring(+ + -)
            if (dif[0] < 0){
                // construct new triangle
                for (int j = 0; j < tri.line_seg.size(); j++) {
                    of_triangle tmpTri(tri.p[0], tri.line_seg[j].p[0], tri.line_seg[j].p[1]);
                    tmpTri.copyNormal(tri.normal);
                    Tlist.push_back(tmpTri);
                    T_num++;
                }
            }
            else if (dif[1] < 0){
                // construct new triangle
                for (int j = 0; j < tri.line_seg.size(); j++) {
                    of_triangle tmpTri(tri.p[1], tri.line_seg[j].p[0], tri.line_seg[j].p[1]);
                    tmpTri.copyNormal(tri.normal);
                    Tlist.push_back(tmpTri);
                    T_num++;
                }
            }
            else if (dif[2] < 0){
                // construct new triangle
                for (int j = 0; j < tri.line_seg.size(); j++) {
                    of_triangle tmpTri(tri.p[2], tri.line_seg[j].p[0], tri.line_seg[j].p[1]);
                    tmpTri.copyNormal(tri.normal);
                    Tlist.push_back(tmpTri);
                    T_num++;
                }
            }
            else
                printf("Err\n");
            
        }
        else if (dif[0] * dif[1] * dif[2] > 0){ // one vertex in ring(- - +)
            
            if (dif[0] > 0){
                if (tri.line_seg[0].p[0].y < tri.line_seg[0].p[1].y) { // convex
                    of_triangle tmpTri2(tri.p[1], tri.p[2], tri.line_seg[0].p[0]);
                    tmpTri2.copyNormal(tri.normal);
                    Tlist.push_back(tmpTri2);
                    T_num++;
                    
                    for (int j = 0; j < tri.line_seg.size(); j++) {
                        of_triangle tmpTri(tri.p[2], tri.line_seg[j].p[0], tri.line_seg[j].p[1]);
                        tmpTri.copyNormal(tri.normal);
                        Tlist.push_back(tmpTri);
                        T_num++;
                    }
                }
                else{
                    // p[0], p[1], line_seg[0].p[0] in the same line
                    if ((tri.p[0] - tri.p[1]).dot(tri.p[0] - tri.line_seg[0].p[0]) == 0) {
                        for (int i = 0; i < tri.line_seg.size(); i++) {
                            bool isMin = false;
                            if (tri.line_seg[i].p[0].y > tri.line_seg[i].p[1].y) {
                                of_triangle tmpTri(tri.p[1], tri.line_seg[i].p[0], tri.line_seg[i].p[1]);
                                tmpTri.copyNormal(tri.normal);
                                Tlist.push_back(tmpTri);
                                T_num++;
                            }
                            else{
                                if (!isMin) {
                                    isMin = !isMin;
                                    of_triangle tmpTri(tri.p[1], tri.p[2], tri.line_seg[i].p[0]);
                                    tmpTri.copyNormal(tri.normal);
                                    Tlist.push_back(tmpTri);
                                    T_num++;
                                }
                                of_triangle tmpTri(tri.p[2], tri.line_seg[i].p[0], tri.line_seg[i].p[1]);
                                tmpTri.copyNormal(tri.normal);
                                Tlist.push_back(tmpTri);
                                T_num++;
                            }
                        }
                    }
                    // p[0], p[2], line_seg[0].p[0] in the same line
                    else{
                        for (int i = 0; i < tri.line_seg.size(); i++) {
                            bool isMin = false;
                            if (tri.line_seg[i].p[0].y > tri.line_seg[i].p[1].y) {
                                of_triangle tmpTri(tri.p[2], tri.line_seg[i].p[0], tri.line_seg[i].p[1]);
                                tmpTri.copyNormal(tri.normal);
                                Tlist.push_back(tmpTri);
                                T_num++;
                            }
                            else{
                                if (!isMin) {
                                    isMin = !isMin;
                                    of_triangle tmpTri(tri.p[2], tri.p[1], tri.line_seg[i].p[0]);
                                    tmpTri.copyNormal(tri.normal);
                                    Tlist.push_back(tmpTri);
                                    T_num++;
                                }
                                of_triangle tmpTri(tri.p[1], tri.line_seg[i].p[0], tri.line_seg[i].p[1]);
                                tmpTri.copyNormal(tri.normal);
                                Tlist.push_back(tmpTri);
                                T_num++;
                            }
                        }
                    }
                }
            }
            else if (dif[1] > 0){
                if (tri.line_seg[0].p[0].y < tri.line_seg[0].p[1].y) { // convex
                    of_triangle tmpTri2(tri.p[2], tri.p[0], tri.line_seg[0].p[0]);
                    tmpTri2.copyNormal(tri.normal);
                    Tlist.push_back(tmpTri2);
                    T_num++;
                    
                    for (int j = 0; j < tri.line_seg.size(); j++) {
                        of_triangle tmpTri(tri.p[0], tri.line_seg[j].p[0], tri.line_seg[j].p[1]);
                        tmpTri.copyNormal(tri.normal);
                        Tlist.push_back(tmpTri);
                        T_num++;
                    }
                }
                else{
                    // p[1], p[2], line_seg[0].p[0] in the same line
                    if ((tri.p[1] - tri.p[2]).dot(tri.p[1] - tri.line_seg[0].p[0]) == 0) {
                        for (int i = 0; i < tri.line_seg.size(); i++) {
                            bool isMin = false;
                            if (tri.line_seg[i].p[0].y > tri.line_seg[i].p[1].y) {
                                of_triangle tmpTri(tri.p[2], tri.line_seg[i].p[0], tri.line_seg[i].p[1]);
                                tmpTri.copyNormal(tri.normal);
                                Tlist.push_back(tmpTri);
                                T_num++;
                            }
                            else{
                                if (!isMin) {
                                    isMin = !isMin;
                                    of_triangle tmpTri(tri.p[2], tri.p[0], tri.line_seg[i].p[0]);
                                    tmpTri.copyNormal(tri.normal);
                                    Tlist.push_back(tmpTri);
                                    T_num++;
                                }
                                of_triangle tmpTri(tri.p[0], tri.line_seg[i].p[0], tri.line_seg[i].p[1]);
                                tmpTri.copyNormal(tri.normal);
                                Tlist.push_back(tmpTri);
                                T_num++;
                            }
                        }
                    }
                    // p[1], p[0], line_seg[0].p[0] in the same line
                    else{
                        for (int i = 0; i < tri.line_seg.size(); i++) {
                            bool isMin = false;
                            if (tri.line_seg[i].p[0].y > tri.line_seg[i].p[1].y) {
                                of_triangle tmpTri(tri.p[0], tri.line_seg[i].p[0], tri.line_seg[i].p[1]);
                                tmpTri.copyNormal(tri.normal);
                                Tlist.push_back(tmpTri);
                                T_num++;
                            }
                            else{
                                if (!isMin) {
                                    isMin = !isMin;
                                    of_triangle tmpTri(tri.p[0], tri.p[2], tri.line_seg[i].p[0]);
                                    tmpTri.copyNormal(tri.normal);
                                    Tlist.push_back(tmpTri);
                                    T_num++;
                                }
                                of_triangle tmpTri(tri.p[2], tri.line_seg[i].p[0], tri.line_seg[i].p[1]);
                                tmpTri.copyNormal(tri.normal);
                                Tlist.push_back(tmpTri);
                                T_num++;
                            }
                        }
                    }
                }
            }
            else if (dif[2] > 0){
                if (tri.line_seg[0].p[0].y < tri.line_seg[0].p[1].y) { // convex
                    of_triangle tmpTri2(tri.p[0], tri.p[1], tri.line_seg[0].p[0]);
                    tmpTri2.copyNormal(tri.normal);
                    Tlist.push_back(tmpTri2);
                    T_num++;
                    
                    for (int j = 0; j < tri.line_seg.size(); j++) {
                        of_triangle tmpTri(tri.p[1], tri.line_seg[j].p[0], tri.line_seg[j].p[1]);
                        tmpTri.copyNormal(tri.normal);
                        Tlist.push_back(tmpTri);
                        T_num++;
                    }
                }
                else{
                    // p[2], p[0], line_seg[0].p[0] in the same line
                    if ((tri.p[2] - tri.p[0]).dot(tri.p[2] - tri.line_seg[0].p[0]) == 0) {
                        for (int i = 0; i < tri.line_seg.size(); i++) {
                            bool isMin = false;
                            if (tri.line_seg[i].p[0].y > tri.line_seg[i].p[1].y) {
                                of_triangle tmpTri(tri.p[0], tri.line_seg[i].p[0], tri.line_seg[i].p[1]);
                                tmpTri.copyNormal(tri.normal);
                                Tlist.push_back(tmpTri);
                                T_num++;
                            }
                            else{
                                if (!isMin) {
                                    isMin = !isMin;
                                    of_triangle tmpTri(tri.p[0], tri.p[1], tri.line_seg[i].p[0]);
                                    tmpTri.copyNormal(tri.normal);
                                    Tlist.push_back(tmpTri);
                                    T_num++;
                                }
                                of_triangle tmpTri(tri.p[1], tri.line_seg[i].p[0], tri.line_seg[i].p[1]);
                                tmpTri.copyNormal(tri.normal);
                                Tlist.push_back(tmpTri);
                                T_num++;
                            }
                        }
                    }
                    // p[2], p[1], line_seg[0].p[0] in the same line
                    else{
                        for (int i = 0; i < tri.line_seg.size(); i++) {
                            bool isMin = false;
                            if (tri.line_seg[i].p[0].y > tri.line_seg[i].p[1].y) {
                                of_triangle tmpTri(tri.p[1], tri.line_seg[i].p[0], tri.line_seg[i].p[1]);
                                tmpTri.copyNormal(tri.normal);
                                Tlist.push_back(tmpTri);
                                T_num++;
                            }
                            else{
                                if (!isMin) {
                                    isMin = !isMin;
                                    of_triangle tmpTri(tri.p[1], tri.p[0], tri.line_seg[i].p[0]);
                                    tmpTri.copyNormal(tri.normal);
                                    Tlist.push_back(tmpTri);
                                    T_num++;
                                }
                                of_triangle tmpTri(tri.p[0], tri.line_seg[i].p[0], tri.line_seg[i].p[1]);
                                tmpTri.copyNormal(tri.normal);
                                Tlist.push_back(tmpTri);
                                T_num++;
                            }
                        }
                    }
                }
                
            }
            else
                printf("Err\n");
            /*
             if (dif[0] > 0){
             if ((tri.p[0] - tri.p[1]).dot(tri.p[0] - tri.line_seg[0].p[0]) == 0) {
             of_triangle tmpTri(tri.p[1], tri.p[2], tri.line_seg[0].p[0]);
             tmpTri.copyNormal(tri.normal);
             Tlist.push_back(tmpTri);
             T_num++;
             
             of_triangle tmpTri2(tri.p[2], tri.line_seg[0].p[0], tri.line_seg[tri.line_seg.size()-1].p[1]);
             tmpTri2.copyNormal(tri.normal);
             Tlist.push_back(tmpTri2);
             T_num++;
             }
             else{
             of_triangle tmpTri(tri.p[1], tri.p[2], tri.line_seg[tri.line_seg.size()-1].p[1]);
             tmpTri.copyNormal(tri.normal);
             Tlist.push_back(tmpTri);
             T_num++;
             
             of_triangle tmpTri2(tri.p[2], tri.line_seg[0].p[0], tri.line_seg[0].p[0]);
             tmpTri2.copyNormal(tri.normal);
             Tlist.push_back(tmpTri2);
             T_num++;
             }
             }
             else if (dif[1] > 0){
             if ((tri.p[1] - tri.p[2]).dot(tri.p[1] - tri.line_seg[0].p[0]) == 0) {
             of_triangle tmpTri(tri.p[2], tri.p[0], tri.line_seg[0].p[0]);
             tmpTri.copyNormal(tri.normal);
             Tlist.push_back(tmpTri);
             T_num++;
             
             of_triangle tmpTri2(tri.p[0], tri.line_seg[0].p[0], tri.line_seg[tri.line_seg.size()-1].p[1]);
             tmpTri2.copyNormal(tri.normal);
             Tlist.push_back(tmpTri2);
             T_num++;
             }
             else{
             of_triangle tmpTri(tri.p[2], tri.p[0], tri.line_seg[tri.line_seg.size()-1].p[1]);
             tmpTri.copyNormal(tri.normal);
             Tlist.push_back(tmpTri);
             T_num++;
             
             of_triangle tmpTri2(tri.p[0], tri.line_seg[0].p[0], tri.line_seg[0].p[0]);
             tmpTri2.copyNormal(tri.normal);
             Tlist.push_back(tmpTri2);
             T_num++;
             }
             }
             else if (dif[2] > 0){
             if ((tri.p[2] - tri.p[0]).dot(tri.p[2] - tri.line_seg[0].p[0]) == 0) {
             of_triangle tmpTri(tri.p[0], tri.p[1], tri.line_seg[0].p[0]);
             tmpTri.copyNormal(tri.normal);
             Tlist.push_back(tmpTri);
             T_num++;
             
             of_triangle tmpTri2(tri.p[1], tri.line_seg[0].p[0], tri.line_seg[tri.line_seg.size()-1].p[1]);
             tmpTri2.copyNormal(tri.normal);
             Tlist.push_back(tmpTri2);
             T_num++;
             }
             else{
             of_triangle tmpTri(tri.p[0], tri.p[1], tri.line_seg[tri.line_seg.size()-1].p[1]);
             tmpTri.copyNormal(tri.normal);
             Tlist.push_back(tmpTri);
             T_num++;
             
             of_triangle tmpTri2(tri.p[1], tri.line_seg[0].p[0], tri.line_seg[0].p[0]);
             tmpTri2.copyNormal(tri.normal);
             Tlist.push_back(tmpTri2);
             T_num++;
             }
             }
             else
             printf("Err\n");
             */
        }
    }
    
#pragma mark TODO
    // calculate extrusing plane
    // Get two end of ring
}

/*

void testApp::removetriBelongToRing(){
    
    float level = Tlist[triBelongToRing[0]].line_seg[0].p[0].y;
    
    // remove the triangle in the ring
    // There would be some bug is the shape is not convex
    for (vector<of_triangle>::iterator it = Tlist.begin(); it != Tlist.end();) {
        if ((*it).p[0].y > level && (*it).p[1].y > level && (*it).p[2].y > level)
            it = Tlist.erase(it);
        else
            it++;
    }
    
    triBelongToRing.clear();
}

*/

void testApp::removeRing(){
    std::sort(triBelongToRing.begin(), triBelongToRing.end());
    triBelongToRing.erase(unique(triBelongToRing.begin(), triBelongToRing.end()), triBelongToRing.end());
    for (int i = triBelongToRing.size()-1; i >= 0; i--) {
        Tlist.erase(Tlist.begin() + triBelongToRing[i]);
    }
}


