#include "testApp.h"

bool release = false;
bool draw_c = false;
bool elevated_T = false;
ofLight l1;

void testApp::elevate()
{
    ofVec3f v;
}


void testApp::add_fan_triangle_1(vector<of_edge>& collect_edge,ofPoint center)
{
    ofPoint tmpP[3];
    tmpP[0] = center;
    for(int i = 0; i<collect_edge.size(); ++i)
    {
        tmpP[1] = collect_edge[i].p[0];
        tmpP[2] = collect_edge[i].p[1];
        Tlist.push_back(tmpP);
        Tlist[T_num].counter[0] = 1;
        Tlist[T_num].type = 3;
        ++T_num;
    }
}

void testApp::add_fan_triangle_2(ofPoint P1,ofPoint P2,ofPoint P3)
{
    ofPoint tmpP[3];
    tmpP[0] = P1;
    tmpP[1] = P2;
    tmpP[2] = P3;
    Tlist.push_back(tmpP);
    Tlist[T_num].counter[0] = 1;
    Tlist[T_num].counter[1] = 1;
    ofPoint c_p[2];
    c_p[0] = P1;
    c_p[1] = P2;
    Tlist[T_num].chordal_axis.push_back(c_p);
    Tlist[T_num].type = 3;
    ++T_num;
}

void testApp::add_fan_triangle_3(ofPoint P1,ofPoint P2,ofPoint P3)
{
    ofPoint tmpP[3];
    tmpP[0] = P1;
    tmpP[1] = P2;
    tmpP[2] = P3;
    Tlist.push_back(tmpP);
    Tlist[T_num].counter[0] = 1;
    Tlist[T_num].type = 3;
    ++T_num;
}

bool testApp::check_same_egde(ofPoint& e1p1,ofPoint& e1p2,ofPoint& e2p1,ofPoint& e2p2)
{
    return ( (e1p1 ==  e2p1 && e1p2 == e2p2) || (e1p1 ==  e2p2 && e1p2 == e2p1) );
}

bool testApp::check_advance(vector<ofPoint>& collect_point,ofPoint& center,float r)
{
    for(int i =0; i < collect_point.size(); ++i)
    {
        if( center.distance(collect_point[i]) > r )
            return false;
    }
    return true;
}

void testApp::prune_1()
{
    ofPoint p_edge[2],center,center2;
    float r;
    for( vector<of_triangle>::iterator i=Tlist.begin(); i!=Tlist.end(); )
    {
        vector<ofPoint> collect_point;
        vector<of_edge> collect_edge;
        if(i->type == 2)
        {
            ofPoint tmp_add_edge[2];
            if(i->counter[0]==2)
            {
                tmp_add_edge[0] = i->p[0];
                tmp_add_edge[1] = i->p[1];
                collect_edge.push_back(tmp_add_edge);
                tmp_add_edge[1] = i->p[2];
                collect_edge.push_back(tmp_add_edge);
                
                collect_point.push_back(i->p[0]);
                p_edge[0] = i->p[1];
                p_edge[1] = i->p[2];
                
            }
            else if(i->counter[1]==2)
            {
                tmp_add_edge[0] = i->p[1];
                tmp_add_edge[1] = i->p[0];
                collect_edge.push_back(tmp_add_edge);
                tmp_add_edge[1] = i->p[2];
                collect_edge.push_back(tmp_add_edge);
                
                collect_point.push_back(i->p[1]);
                p_edge[0] = i->p[0];
                p_edge[1] = i->p[2];

            }
            else
            {
                tmp_add_edge[0] = i->p[2];
                tmp_add_edge[1] = i->p[0];
                collect_edge.push_back(tmp_add_edge);
                tmp_add_edge[1] = i->p[1];
                collect_edge.push_back(tmp_add_edge);
                
                collect_point.push_back(i->p[2]);
                p_edge[0] = i->p[0];
                p_edge[1] = i->p[1];
            }
            
            r = ( p_edge[0].distance(p_edge[1]) )/2.0f;
            center = (p_edge[0] + p_edge[1])/2;
            i = Tlist.erase(i);
            --T_num;
            bool advanced = true;//check_advance(collect_point, center, r);
            
            while(advanced)
            {
                for(vector<of_triangle>::iterator j=Tlist.begin(); j!=Tlist.end();)
                {
                    if(
                       check_same_egde(p_edge[0], p_edge[1], j->p[0], j->p[1])
                       )
                    {
                        if(j->type == 0)
                        {
                            
                            center = ( j->p[0] + j->p[1] + j->p[2] )/3;
                            ofPoint find_c  = (j->p[0] + j->p[1])/2 ;
                            for(vector<of_edge>::iterator k=j->chordal_axis.begin(); k!=j->chordal_axis.end();)
                            {
                                if(check_same_egde(center, find_c, k->p[0], k->p[1]))
                                    k = j->chordal_axis.erase(k);
                                else
                                    ++k;
                            }
                            
                            add_fan_triangle_1(collect_edge,center);
                             /*
                             add_fan_triangle_2(center, ( j->p[1] + j->p[2] )/2,j->p[1]);
                             add_fan_triangle_2(center, ( j->p[1] + j->p[2] )/2,j->p[2]);
                             
                             add_fan_triangle_2(center, ( j->p[2] + j->p[0] )/2,j->p[2]);
                             add_fan_triangle_2(center, ( j->p[2] + j->p[0] )/2,j->p[0]);
                             */
                              
                            /*if(i == j)
                                ++i;
                            Tlist.erase(j);
                            --T_num;*/
                            advanced = false;
                            break;
                        }
                        if(j->type == 1)
                        {
                            p_edge[0] = j->p[2];
                            if(j->counter[0] == 1)
                            {
                                p_edge[1] = j->p[0];
                                collect_edge.push_back(p_edge);
                                p_edge[1] = j->p[1];
                                collect_point.push_back(j->p[0]);
                            }
                            else
                            {
                                p_edge[1] = j->p[1];
                                collect_edge.push_back(p_edge);
                                p_edge[1] = j->p[0];
                                collect_point.push_back(j->p[1]);
                            }
                            
                            center = ( p_edge[0] + p_edge[1] )/2;
                            r = ( p_edge[0].distance(p_edge[1]) )/2;
                            advanced = check_advance(collect_point, center, r);
                            
                            
                            if(!advanced)
                                add_fan_triangle_1(collect_edge, center);
                            
                            if(i == j)
                                ++i;
                            Tlist.erase(j);
                            --T_num;
                            break;
                        }
                    }
                    else if(
                       check_same_egde(p_edge[0], p_edge[1], j->p[1], j->p[2])
                       )
                    {
                        if(j->type == 0)
                        {
                            center = ( j->p[0] + j->p[1] + j->p[2] )/3;
                            ofPoint find_c  = (j->p[1] + j->p[2])/2 ;
                            for(vector<of_edge>::iterator k=j->chordal_axis.begin(); k!=j->chordal_axis.end();)
                            {
                                if(check_same_egde(center, find_c, k->p[0], k->p[1]))
                                    k = j->chordal_axis.erase(k);
                                else
                                    ++k;
                            }
                            
                             add_fan_triangle_1(collect_edge,center);
                             /*
                             add_fan_triangle_2(center, ( j->p[0] + j->p[1] )/2,j->p[0]);
                             add_fan_triangle_2(center, ( j->p[0] + j->p[1] )/2,j->p[1]);
                             
                             add_fan_triangle_2(center, ( j->p[2] + j->p[0] )/2,j->p[2]);
                             add_fan_triangle_2(center, ( j->p[2] + j->p[0] )/2,j->p[0]);
                             */
                            
                            /*if(i == j)
                                ++i;
                            Tlist.erase(j);
                            --T_num;*/
                            advanced = false;
                            break;
                        }
                        if(j->type == 1)
                        {
                            p_edge[0] = j->p[0];
                            if(j->counter[1] == 1)
                            {
                                p_edge[1] = j->p[1];
                                collect_edge.push_back(p_edge);
                                p_edge[1] = j->p[2];
                                collect_point.push_back(j->p[1]);
                            }
                            else
                            {
                                p_edge[1] = j->p[2];
                                collect_edge.push_back(p_edge);
                                p_edge[1] = j->p[1];
                                collect_point.push_back(j->p[2]);
                            }
                            
                            center = ( p_edge[0] + p_edge[1] )/2;
                            r = ( p_edge[0].distance(p_edge[1]) )/2;
                            advanced = check_advance(collect_point, center, r);
                            
                            if(!advanced)
                                add_fan_triangle_1(collect_edge, center);
                            
                            if(i == j)
                                ++i;
                            Tlist.erase(j);
                            --T_num;
                            break;
                        }
                    }
                    else if(
                       check_same_egde(p_edge[0], p_edge[1], j->p[2], j->p[0])
                       )
                    {
                        if(j->type == 0)
                        {
                            
                             center = ( j->p[0] + j->p[1] + j->p[2] )/3;
                            ofPoint find_c  = (j->p[2] + j->p[0])/2 ;
                            for(vector<of_edge>::iterator k=j->chordal_axis.begin(); k!=j->chordal_axis.end();)
                            {
                                if(check_same_egde(center, find_c, k->p[0], k->p[1]))
                                    k = j->chordal_axis.erase(k);
                                else
                                    ++k;
                            }
                             
                             add_fan_triangle_1(collect_edge,center);
                             /*
                             add_fan_triangle_2(center, ( j->p[0] + j->p[1] )/2,j->p[0]);
                             add_fan_triangle_2(center, ( j->p[0] + j->p[1] )/2,j->p[1]);
                             
                             add_fan_triangle_2(center, ( j->p[1] + j->p[2] )/2,j->p[1]);
                             add_fan_triangle_2(center, ( j->p[1] + j->p[2] )/2,j->p[2]);
                            */
                            
                            /*if(i == j)
                                ++i;
                            Tlist.erase(j);
                            --T_num;*/
                            advanced = false;
                            break;
                        }
                        if(j->type == 1)
                        {
                            p_edge[0] = j->p[1];
                            if(j->counter[2] == 1)
                            {
                                p_edge[1] = j->p[2];
                                collect_edge.push_back(p_edge);
                                p_edge[1] = j->p[0];
                                collect_point.push_back(j->p[2]);
                            }
                            else
                            {
                                p_edge[1] = j->p[0];
                                collect_edge.push_back(p_edge);
                                p_edge[1] = j->p[2];
                                collect_point.push_back(j->p[0]);
                            }
                            
                            center = ( p_edge[0] + p_edge[1] )/2;
                            r = ( p_edge[0].distance(p_edge[1]) )/2;
                            advanced = check_advance(collect_point, center, r);
                            
                            if(!advanced)
                                add_fan_triangle_1(collect_edge, center);
                            
                            if(i == j)
                                ++i;
                            Tlist.erase(j);
                            --T_num;
                            break;
                        }
                        

                    }
                    
                    else
                        ++j;
                    
                    if(j == Tlist.end())
                        advanced = false;
                }
                
            }
            
        }
        else
            ++i;
    }
}

void testApp::prune_2()
{
    ofPoint center,center2;
    for( int i=0 ; i<T_num ; ++i )
    {
        if( Tlist[i].type == 1 )
        {
            if(Tlist[i].counter[0] == 0)
            {
                center = ( Tlist[i].p[0] + Tlist[i].p[1] )/2;
                center2 =( Tlist[i].p[0] + Tlist[i].p[2] )/2;
                
                add_fan_triangle_2(center, center2, Tlist[i].p[0]);
                add_fan_triangle_2(center, center2, Tlist[i].p[1]);
                
                add_fan_triangle_3(center2, Tlist[i].p[1] , Tlist[i].p[2] );
            }
            
            if(Tlist[i].counter[1] == 0)
            {
                center = ( Tlist[i].p[1] + Tlist[i].p[0] )/2;
                center2 =( Tlist[i].p[1] + Tlist[i].p[2] )/2;
                
                add_fan_triangle_2(center, center2, Tlist[i].p[1]);
                add_fan_triangle_2(center, center2, Tlist[i].p[0]);
                
                add_fan_triangle_3(center2, Tlist[i].p[0] , Tlist[i].p[2] );
            }
            
            if(Tlist[i].counter[2] == 0)
            {
                center = ( Tlist[i].p[2] + Tlist[i].p[0] )/2;
                center2 =( Tlist[i].p[2] + Tlist[i].p[1] )/2;
                
                add_fan_triangle_2(center, center2, Tlist[i].p[2]);
                add_fan_triangle_2(center, center2, Tlist[i].p[0]);
                
                add_fan_triangle_3(center2, Tlist[i].p[0] , Tlist[i].p[1] );
            }
        }
        
        if(Tlist[i].type == 0)
        {
            for(int j = 0; j<Tlist[i].chordal_axis.size()  ; ++j)
            {
                ofPoint tmpP[2];
                tmpP[0] = (Tlist[i].p[0]+Tlist[i].p[1]+Tlist[i].p[2])/3;
                
                tmpP[1] = (Tlist[i].p[0]+Tlist[i].p[1])/2;
                if(check_same_egde(tmpP[0],tmpP[1], Tlist[i].chordal_axis[j].p[0], Tlist[i].chordal_axis[j].p[1]))
                {
                    add_fan_triangle_2(tmpP[0], tmpP[1], Tlist[i].p[0]);
                    add_fan_triangle_2(tmpP[0], tmpP[1], Tlist[i].p[1]);
                }
                
                tmpP[1] = (Tlist[i].p[1]+Tlist[i].p[2])/2;
                if(check_same_egde(tmpP[0],tmpP[1], Tlist[i].chordal_axis[j].p[0], Tlist[i].chordal_axis[j].p[1]))
                {
                    add_fan_triangle_2(tmpP[0], tmpP[1], Tlist[i].p[1]);
                    add_fan_triangle_2(tmpP[0], tmpP[1], Tlist[i].p[2]);
                }
                
                tmpP[1] = (Tlist[i].p[2]+Tlist[i].p[0])/2;
                if(check_same_egde(tmpP[0],tmpP[1], Tlist[i].chordal_axis[j].p[0], Tlist[i].chordal_axis[j].p[1]))
                {
                    add_fan_triangle_2(tmpP[0], tmpP[1], Tlist[i].p[2]);
                    add_fan_triangle_2(tmpP[0], tmpP[1], Tlist[i].p[0]);
                }
            }
        }
    }
    
    
    for( vector<of_triangle>::iterator i=Tlist.begin(); i!=Tlist.end(); )
    {
        if(i->type != 3)
        {
            i = Tlist.erase(i);
        }
        else
            ++i;
    }
}




//--------------------------------------------------------------
void testApp::setup(){
    
    ofFloatColor c_l1;
    c_l1.set(20, 20, 20);
    l1.setAmbientColor(c_l1);
    c_l1.set(0, 0, 100);
    l1.setDiffuseColor(c_l1);
    c_l1.set(0, 0, 255);
    l1.setSpecularColor(c_l1);
    l1.setPointLight();
    l1.setPosition(0, 0, -50);
    l1.disable();
    
}

//--------------------------------------------------------------
void testApp::update(){
    
    if(elevated_T)
        l1.enable();
    else
        l1.disable();
    
    
}

//--------------------------------------------------------------
void testApp::draw(){
    
    
    //ofBackgroundGradient(ofColor::white, ofColor(200,200,200), OF_GRADIENT_LINEAR);
    if(elevated_T)
    {
        for(int i =0; i < T_num; ++i)
        {
            Tlist[i].draw_triangle();
        }
    }
    else if(release)
    {
        for(int i =0; i < T_num; ++i)
        {
            ofSetColor(80*(2-Tlist[i].type),80*Tlist[i].type,80*(2-Tlist[i].type));
            Tlist[i].draw_triangle();
            ofSetColor(255, 255, 255);
            Tlist[i].draw_wireframe();
            ofSetColor(0, 0, 0);
            if(draw_c)
                for(int j = 0; j<Tlist[i].chordal_axis.size() ; ++j)
                    Tlist[i].chordal_axis[j].draw();
            //printf("type = %d\n",Tlist[i].type);
        }
        //mesh.draw();
    }
    else
        line.draw();
}

//--------------------------------------------------------------
void testApp::keyPressed(int key){
    
    if(key == 'c')
        draw_c = !draw_c;
    if(key == 'p')
    {
        prune_1();
        prune_2();
    }
    if(key == 'e')
    {
        elevate();
        elevated_T = true;
    }
    
}

//--------------------------------------------------------------
void testApp::keyReleased(int key){
    
}

//--------------------------------------------------------------
void testApp::mouseMoved(int x, int y ){
    
}

//--------------------------------------------------------------
void testApp::mouseDragged(int x, int y, int button){
    
    
    line.addVertex(ofPoint(x, y));
}

//--------------------------------------------------------------
void testApp::mousePressed(int x, int y, int button){
    line.clear();
    Tlist.clear();
    T_num = 0;
    release = false;
    line.addVertex(ofPoint(x, y));
    
}

//--------------------------------------------------------------
void testApp::mouseReleased(int x, int y, int button){
    
    if (line.size() > 2){
        
        release = true;
        
        ofPolyline lineRespaced = line;
        
        // add the last point (so when we resample, it's a closed polygon)
        lineRespaced.addVertex(lineRespaced[0]);
        // resample
        lineRespaced = lineRespaced.getResampledBySpacing(20);
        // I want to make sure the first point and the last point are not the same, since triangle is unhappy:
        lineRespaced.getVertices().erase(lineRespaced.getVertices().begin());
        // if we have a proper set of points, mesh them:
        
        if (lineRespaced.size() > 5){
            // angle constraint = 28
            // size constraint = -1 (don't constraint triangles by size);
            mesh.triangulate(lineRespaced, 0, -1);
            // this is an alternative, constrain on size not angle:
            //mesh.triangulate(lineRespaced, -1, 200);
            // see ofxTriangleMesh.h for info.
            for(int i = 0; i < mesh.nTriangles; ++i)
            {
                ofPoint tmpP[3];
                tmpP[0] = mesh.outputPts[mesh.triangles[i].index[0]];
                tmpP[1] = mesh.outputPts[mesh.triangles[i].index[1]];
                tmpP[2] = mesh.outputPts[mesh.triangles[i].index[2]];
                Tlist.push_back(tmpP);
                int internal_edges = 0;
                //int counter[3]={0,0,0};
                for(int j = 0,k; j< lineRespaced.size() ; ++j)
                {
                    if(j+1 == lineRespaced.size() )
                        k = 0;
                    else
                        k = j+1;
                    //printf("%d %f %f\n",j,lineRespaced[j].x,lineRespaced[j].y);
                    if(
                       check_same_egde(tmpP[0], tmpP[1], lineRespaced[j], lineRespaced[k])
                       )
                        ++ internal_edges,++Tlist[i].counter[0],++Tlist[i].counter[1];
                    if(
                       check_same_egde(tmpP[1], tmpP[2], lineRespaced[j], lineRespaced[k])
                       )
                        ++ internal_edges,++Tlist[i].counter[1],++Tlist[i].counter[2];
                    if(
                       check_same_egde(tmpP[2], tmpP[0], lineRespaced[j], lineRespaced[k])
                       )
                        ++ internal_edges,++Tlist[i].counter[2],++Tlist[i].counter[0];
                    if(internal_edges == 2)
                        break;
                    
                }
                if(internal_edges == 2)
                {
                    int i1,i2,c1;
                    ofPoint tmpE[2];
                    if(Tlist[i].counter[0]==2)
                    {c1=0;i1=1;i2=2;}
                    else if(Tlist[i].counter[1]==2)
                    {c1=1;i1=0;i2=2;}
                    else
                    {c1=2;i1=0;i2=1;}
                    tmpE[0] = tmpP[c1];tmpE[1] = (tmpP[i1]+tmpP[i2])/2;
                    Tlist[T_num].chordal_axis.push_back(tmpE);
                }
                else if(internal_edges == 1)
                {
                    int x1,i2,i1;
                    ofPoint tmpE[2];
                    if(Tlist[i].counter[0]==0)
                    {x1=0;i1=1;i2=2;}
                    else if(Tlist[i].counter[1]==0)
                    {x1=1;i1=0;i2=2;}
                    else
                    {x1=2;i1=0;i2=1;}
                    tmpE[0] = (tmpP[x1]+tmpP[i1])/2;
                    tmpE[1] = (tmpP[x1]+tmpP[i2])/2;
                    Tlist[T_num].chordal_axis.push_back(tmpE);
                }
                else
                {
                    ofPoint tmpE[2];
                    tmpE[0] = (tmpP[0]+tmpP[1]+tmpP[2])/3;
                    tmpE[1] = (tmpP[0]+tmpP[1])/2;
                    Tlist[T_num].chordal_axis.push_back(tmpE);
                    tmpE[1] = (tmpP[1]+tmpP[2])/2;
                    Tlist[T_num].chordal_axis.push_back(tmpE);
                    tmpE[1] = (tmpP[2]+tmpP[0])/2;
                    Tlist[T_num].chordal_axis.push_back(tmpE);
                }
                Tlist[T_num++].type = internal_edges;
                
            }
        }
    }
}

//--------------------------------------------------------------
void testApp::windowResized(int w, int h){
    
}

//--------------------------------------------------------------
void testApp::gotMessage(ofMessage msg){
    
}

//--------------------------------------------------------------
void testApp::dragEvent(ofDragInfo dragInfo){
    
}