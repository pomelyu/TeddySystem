#include "testApp.h"

ofCamera tmp_cam;


void testApp::cast()
{
    for( vector<of_triangle>::iterator i=Tlist.begin(); i!=Tlist.end(); )
    {
        ofVec3f v1 = i->p[2] - i->p[0];
        ofVec3f v2 = i->p[1] - i->p[0];
        v1 = v2.cross(v1);
        if(v1.length() < 0.1)
        {
            i = Tlist.erase(i);
            --T_num;
        }
        else
            ++i;
    }
}

bool testApp::check_notinside(vector<ofPoint>& added_point,ofPoint check_point)
{
    for(int i = 0; i<added_point.size() ; ++i)
    {
        if(check_point == added_point[i])
            return false;
    }
    return true;
}

void testApp::elevate()
{
    for(int i =0; i<Tlist.size(); ++i)
    {
        for(int j = 0; j<3; ++j)
        {
            if(Tlist[i].counter[j] && (Tlist[i].p[j].z == 0))
            {
                ofPoint cur_P = Tlist[i].p[j];
                vector<ofPoint> added_point;
                vector<int> which_T,which_P;
                float total_d=0.0f;
                int total_P=0;
                for(int k = 0; k<Tlist.size(); ++k)
                {
                    if(Tlist[k].p[0] == cur_P)
                    {
                        which_T.push_back(k);
                        which_P.push_back(0);
                        if(!Tlist[k].counter[1] && check_notinside(added_point, Tlist[k].p[1]))
                        {
                            total_d += cur_P.distance(Tlist[k].p[1]);
                            ++total_P;
                            added_point.push_back(Tlist[k].p[1]);
                        }
                        if(!Tlist[k].counter[2] && check_notinside(added_point, Tlist[k].p[2]))
                        {
                            total_d += cur_P.distance(Tlist[k].p[2]);
                            ++total_P;
                            added_point.push_back(Tlist[k].p[2]);
                        }
                        
                    }
                    
                    if(Tlist[k].p[1] == cur_P)
                    {
                        which_T.push_back(k);
                        which_P.push_back(1);
                        if(!Tlist[k].counter[0] && check_notinside(added_point, Tlist[k].p[0]))
                        {
                            total_d += cur_P.distance(Tlist[k].p[0]);
                            ++total_P;
                            added_point.push_back(Tlist[k].p[0]);
                        }
                        if(!Tlist[k].counter[2] && check_notinside(added_point, Tlist[k].p[2]))
                        {
                            total_d += cur_P.distance(Tlist[k].p[2]);
                            ++total_P;
                            added_point.push_back(Tlist[k].p[2]);
                        }
                        

                    }
                    
                    if(Tlist[k].p[2] == cur_P)
                    {
                        which_T.push_back(k);
                        which_P.push_back(2);
                        if(!Tlist[k].counter[0] && check_notinside(added_point, Tlist[k].p[0]))
                        {
                            total_d += cur_P.distance(Tlist[k].p[0]);
                            ++total_P;
                            added_point.push_back(Tlist[k].p[0]);
                        }
                        if(!Tlist[k].counter[1] && check_notinside(added_point, Tlist[k].p[1]))
                        {
                            total_d += cur_P.distance(Tlist[k].p[1]);
                            ++total_P;
                            added_point.push_back(Tlist[k].p[1]);
                        }
                        
                    }
                }
                float shift = (total_P==0)? 0 :total_d/(float)total_P;
                
                ofVec3f normal = ofVec3f(0.0f,0.0f,0.0f);
                
                for(int k = 0; k<which_T.size(); ++k)
                {
                    Tlist[which_T[k]].p[which_P[k]].z = shift;
                    
                    ofVec3f tmp_n1 = Tlist[which_T[k]].p[2] - Tlist[which_T[k]].p[0];
                    ofVec3f tmp_n2 = Tlist[which_T[k]].p[1] - Tlist[which_T[k]].p[0];
                    tmp_n1 = tmp_n2.cross(tmp_n1);
                    if(tmp_n1.z < 0)
                        tmp_n1 = -tmp_n1;
                    tmp_n1 = tmp_n1.normalized();
                    normal += tmp_n1;
                }
                if(normal.z <0)
                    normal = -normal;
                for(int k = 0; k<which_T.size(); ++k)
                {
                    Tlist[which_T[k]].normal[which_P[k]] = normal;
                }
            }
        }
    }
    for(int i =0; i<Tlist.size(); ++i)
    {
        ofVec3f normal = Tlist[i].p[2]- Tlist[i].p[0];
        normal.cross(Tlist[i].p[1]- Tlist[i].p[0]);
        normal = normal.normalized();
        if(normal.z<0)
            normal = -normal;
        for(int j = 0; j<3; ++j)
        {
            if(!Tlist[i].counter[j])
            {
                Tlist[i].normal[j] = normal;
            }
        }
    }
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
            --T_num;
        }
        else
            ++i;
    }
}




//--------------------------------------------------------------
void testApp::setup(){
    
    ofSetSmoothLighting(true);
    pointLight.setDiffuseColor( ofFloatColor(1.f, 0.f, 0.0f) );
    pointLight.setSpecularColor( ofFloatColor(1.0f, 0.f, 0.f));
    pointLight.setPosition(512, 384, 600);
    
    // shininess is a value between 0 - 128, 128 being the most shiny //
	material.setShininess( 50 );
    // the light highlight of the material //
	material.setSpecularColor(ofColor(255, 255, 255, 255));
    
    cam.resetTransform();
    cam.setFov(90);
    cam.clearParent();
    cam.setPosition(512, 384, 400);
    cam.lookAt(ofVec3f(512,384,0));
    //cam.enableOrtho();
    
    tmp_cam.resetTransform();
    tmp_cam.setFov(120);
    tmp_cam.clearParent();
    tmp_cam.setPosition(0, 0, 0);
    tmp_cam.lookAt(ofVec3f(512,384,0));
    
}

//--------------------------------------------------------------
void testApp::update(){
}

//--------------------------------------------------------------
void testApp::draw(){
    
    
    //ofBackgroundGradient(ofColor::white, ofColor(200,200,200), OF_GRADIENT_LINEAR);
    if(elevated_T)
    {
        ofEnableDepthTest();
        
        ofEnableLighting();
        pointLight.enable();
        
        cam.begin();
        material.begin();
        ofSetColor(255,255,255);
        
        for(int i =0; i < T_num; ++i)
        {
            glBegin(GL_TRIANGLES);
                glNormal3f(Tlist[i].normal[0].x,Tlist[i].normal[0].y,Tlist[i].normal[0].z);
                glVertex3f(Tlist[i].p[0].x,Tlist[i].p[0].y,Tlist[i].p[0].z);
                glNormal3f(Tlist[i].normal[1].x,Tlist[i].normal[1].y,Tlist[i].normal[1].z);
                glVertex3f(Tlist[i].p[1].x,Tlist[i].p[1].y,Tlist[i].p[1].z);
                glNormal3f(Tlist[i].normal[2].x,Tlist[i].normal[2].y,Tlist[i].normal[2].z);
                glVertex3f(Tlist[i].p[2].x,Tlist[i].p[2].y,Tlist[i].p[2].z);
            glEnd();
        }
        
        material.end();
        cam.end();
        
        ofDisableLighting();
        ofDisableDepthTest();
        
        //ofSetColor(pointLight.getDiffuseColor());
        //pointLight.draw();
        
        //for(int i =0; i < T_num; ++i)
        //{
            //ofSetColor(80*(2-Tlist[i].type),80*Tlist[i].type,80*(2-Tlist[i].type));
            //Tlist[i].draw_triangle();
            //ofSetColor(255, 255, 255);
            //Tlist[i].draw_wireframe();
        //}
    }
    else if(release)
    {
        //tmp_cam.begin();
        cam.begin();
        for(int i =0; i < T_num; ++i)
        {
            ofSetColor(80*(2-Tlist[i].type),80*Tlist[i].type,128);
            Tlist[i].draw_triangle();
            ofSetColor(255, 255, 255);
            Tlist[i].draw_wireframe();
            ofSetColor(0, 0, 0);
            if(draw_c)
                for(int j = 0; j<Tlist[i].chordal_axis.size() ; ++j)
                    Tlist[i].chordal_axis[j].draw();
            //printf("type = %d\n",Tlist[i].type);
        }
        //cam.draw();
        //mesh.draw();
        //pointLight.draw();
        cam.end();
        //tmp_cam.end();
    }
    else
    {
        //tmp_cam.begin();
        cam.begin();
        line.draw();
        //cam.draw();
        //pointLight.draw();
        cam.end();
        //tmp_cam.end();
    }
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
    if(key == '1')
    {
        for(int i = 0; i<T_num; ++i)
        {
            Tlist[i].p[0] = Tlist[i].p[0].rotate(1, ofVec3f(512,384,0), ofVec3f(0,0,1));
            Tlist[i].p[1] = Tlist[i].p[1].rotate(1, ofVec3f(512,384,0), ofVec3f(0,0,1));
            Tlist[i].p[2] = Tlist[i].p[2].rotate(1, ofVec3f(512,384,0), ofVec3f(0,0,1));
            
            Tlist[i].normal[0] = Tlist[i].normal[0].rotate(-1, ofVec3f(0,0,1));
            Tlist[i].normal[1] = Tlist[i].normal[1].rotate(-1, ofVec3f(0,0,1));
            Tlist[i].normal[2] = Tlist[i].normal[2].rotate(-1, ofVec3f(0,0,1));
                
        }
    }
    if(key == '2')
    {
        for(int i = 0; i<T_num; ++i)
        {
            Tlist[i].p[0] = Tlist[i].p[0].rotate(1, ofVec3f(512,384,0), ofVec3f(0,1,0));
            Tlist[i].p[1] = Tlist[i].p[1].rotate(1, ofVec3f(512,384,0), ofVec3f(0,1,0));
            Tlist[i].p[2] = Tlist[i].p[2].rotate(1, ofVec3f(512,384,0), ofVec3f(0,1,0));
            
            Tlist[i].normal[0] = Tlist[i].normal[0].rotate(-1, ofVec3f(0,1,0));
            Tlist[i].normal[1] = Tlist[i].normal[1].rotate(-1, ofVec3f(0,1,0));
            Tlist[i].normal[2] = Tlist[i].normal[2].rotate(-1, ofVec3f(0,1,0));
            
        }
    }
    if(key == '3')
    {
        for(int i = 0; i<T_num; ++i)
        {
            Tlist[i].p[0] = Tlist[i].p[0].rotate(1, ofVec3f(512,384,0), ofVec3f(1,0,0));
            Tlist[i].p[1] = Tlist[i].p[1].rotate(1, ofVec3f(512,384,0), ofVec3f(1,0,0));
            Tlist[i].p[2] = Tlist[i].p[2].rotate(1, ofVec3f(512,384,0), ofVec3f(1,0,0));
            
            Tlist[i].normal[0] = Tlist[i].normal[0].rotate(-1, ofVec3f(1,0,0));
            Tlist[i].normal[1] = Tlist[i].normal[1].rotate(-1, ofVec3f(1,0,0));
            Tlist[i].normal[2] = Tlist[i].normal[2].rotate(-1, ofVec3f(1,0,0));
            
        }
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
    
    
    line.addVertex(ofPoint(x, 768-y));
}

//--------------------------------------------------------------
void testApp::mousePressed(int x, int y, int button){
    line.clear();
    Tlist.clear();
    T_num = 0;
    release = false;
    elevated_T = false;
    line.addVertex(ofPoint(x, 768-y));
    
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
            
            cast();
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