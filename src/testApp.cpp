#include "testApp.h"

ofCamera tmp_cam;
float angle = 5.0f;

void add_cut_point(vector<ofPoint>& cut_point,ofPoint point)
{
    for (int i=0; i<cut_point.size(); ++i) {
        if(cut_point[i] == point)
            return;
    }
    cut_point.push_back(point);
}

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

vector<of_edge> cut_line;
vector<ofPoint> cut_point;
int index_i;

void testApp::cut()
{
    vector<int> which_T;
    vector<of_triangle> tmp_t;
    for(int j =0; j<Tlist.size(); ++j)
    {
        which_T.push_back(1);
    }
    for(int i =0; i<cut_point.size(); ++i)
    //int i=index_i;
    //if(i<cut_point.size())
    {
        
        ofPoint tmp_p = cut_point[i];
        for(int j =0; j<Tlist.size(); ++j)
        {
            if(which_T[j]==1)
            {
                if(Tlist[j].p[0] == tmp_p)
                {
                    add_cut_point(cut_point, Tlist[j].p[1]);
                    add_cut_point(cut_point, Tlist[j].p[2]);
                    which_T[j] = 0;
                }
                else if (Tlist[j].p[1] == tmp_p)
                {
                    add_cut_point(cut_point, Tlist[j].p[0]);
                    add_cut_point(cut_point, Tlist[j].p[2]);
                    which_T[j] = 0;
                }
                else if (Tlist[j].p[2] == tmp_p)
                {
                    add_cut_point(cut_point, Tlist[j].p[1]);
                    add_cut_point(cut_point, Tlist[j].p[1]);
                    which_T[j] = 0;

                }
            }
        }
        //++index_i;
    }
    for(int j =0; j<Tlist.size(); ++j)
    {
        if(which_T[j])
        {
            tmp_t.push_back(Tlist[j]);
        }
    }
    Tlist = tmp_t;
    T_num = tmp_t.size();
}


/*void testApp::cut()
{
    for (vector<ofPoint>::iterator i = cut_point.begin(); cut_point.size()!=0; )
    //if(cut_point.size() != 0)
    {
        //vector<ofPoint>::iterator i = cut_point.begin();
        
        ofPoint tmp_p = *i;
        cut_point.erase(i);
        for(vector<of_triangle>::iterator j = Tlist.begin(); j!=Tlist.end();)
        {
            if(j->p[0] == tmp_p)
            {
                add_cut_point(cut_point, j->p[1]);
                add_cut_point(cut_point, j->p[2]);
                j = Tlist.erase(j);
            }
            else if (j->p[1] == tmp_p)
            {
                add_cut_point(cut_point, j->p[0]);
                add_cut_point(cut_point, j->p[2]);
                j = Tlist.erase(j);
            }
            else if (j->p[2] == tmp_p)
            {
                add_cut_point(cut_point, j->p[0]);
                add_cut_point(cut_point, j->p[1]);
                j = Tlist.erase(j);
            }
            else
                ++j;
        }
        i = cut_point.begin();
    }
}
*/

bool check_point_same(ofPoint& p1,ofPoint& p2)
{
    float dx = p1.x-p2.x;
    float dy = p1.y-p2.y;
    float dz = p1.z-p2.z;
    if(abs(dx) < 0.01 && abs(dy) < 0.01 && abs(dz) < 0.01)
    {
        return true;
    }
    return false;
}
void testApp::cut_construct(ofVec3f& plane_normal)
{
    //int added_edge;
    /*vector<int> used;
    for(int j =0; j<cut_line.size();++j)
    {
        used.push_back(0);
    }
    ofPolyline cut_line_construct;
    cut_line_construct.addVertex(cut_line[0].p[0]);
    for(int i =0; i!= cut_line.size();)
    {
        bool find = false;
        for(int j =0; j<cut_line.size();++j)
        {
            if(check_point_same( cut_line[j].p[0],cut_line_construct[i]) && used[j]==0 )
            {
                cut_line_construct.addVertex(cut_line[j].p[1]);
                find = true;
                used[j] = j+1;
                ++i;
                //break;
            }
            else if (check_point_same( cut_line[j].p[1],cut_line_construct[i]) && used[j]==0 )
            {
                cut_line_construct.addVertex(cut_line[j].p[0]);
                used[j] = j+1;
                find = true;
                ++i;
                //break;
            }
        }
        if(!find)
            ++i;
    }
    ofxTriangleMesh mesh_construct;
    mesh_construct.triangulate(cut_line_construct, 0, -1);
    
    for(int i = 0; i < mesh_construct.nTriangles; ++i)
    {
        ofPoint tmpP[3];
        tmpP[0] = mesh_construct.outputPts[mesh_construct.triangles[i].index[0]];
        tmpP[1] = mesh_construct.outputPts[mesh_construct.triangles[i].index[1]];
        tmpP[2] = mesh_construct.outputPts[mesh_construct.triangles[i].index[2]];
        Tlist.push_back(tmpP);
        
        ofVec3f v1= tmpP[2] - tmpP[0];
        ofVec3f v2= tmpP[2] - tmpP[0];
        v2.cross(v1);
        v2.normalize();
        if(v2.angle(plane_normal) < 90)
            v2=-v2;
        v2 = ofVec3f(-1,0,0);
        Tlist.back().normal[0] = v2;
        Tlist.back().normal[1] = v2;
        Tlist.back().normal[2] = v2;
    }*/
    
    float cx=0.0f,cy=0.0f,cz=0.0f;
    for(int j =0; j<cut_line.size();++j)
    {
        cx+=cut_line[j].p[0].x;
        cx+=cut_line[j].p[1].x;
        cy+=cut_line[j].p[0].y;
        cy+=cut_line[j].p[1].y;
        cz+=cut_line[j].p[0].z;
        cz+=cut_line[j].p[1].z;
    }
    cx = cx/cut_line.size()/2.0f;
    cy = cy/cut_line.size()/2.0f;
    cz = cz/cut_line.size()/2.0f;
    ofPoint cut_center = ofPoint(cx,cy,cz);
    
    for(int i = 0; i < cut_line.size(); ++i)
    {
        ofPoint tmpP[3];
        tmpP[0] = cut_line[i].p[0];
        tmpP[1] = cut_line[i].p[1];
        tmpP[2] = cut_center;
        Tlist.push_back(tmpP);
        
        ofVec3f v1= tmpP[2] - tmpP[0];
        ofVec3f v2= tmpP[2] - tmpP[0];
        v2.cross(v1);
        v2.normalize();
        if(v2.angle(plane_normal) < 90)
            v2=-v2;
        v2 = ofVec3f(-1,0,0);
        Tlist.back().normal[0] = v2;
        Tlist.back().normal[1] = v2;
        Tlist.back().normal[2] = v2;
    }
    
    
    T_num =  Tlist.size();
    cast();
}

void testApp::cut_plane()
{
    cut_point.clear();
    cut_line.clear();
    index_i = 0;
    
    vector<of_triangle> more_to_add;

    
    vector<int> which_intersect;
    for(int i = 0; i<Tlist.size(); ++i)
    {
        which_intersect.push_back(1);
    }
    
    //ofPolyline lineRespaced = line;
    ofPolyline lineRespaced ;
    //vector<ofPoint> cut_point;
    //lineRespaced = lineRespaced.getResampledBySpacing(30);
    ofPoint original_cam = ofPoint(512,384,400);
    //for(int i =1; i <lineRespaced.size(); ++i)
    //{
    lineRespaced.addVertex(line[0]);
    lineRespaced.addVertex(line[line.size()-1]);
    int i =1;
            ofVec3f v1 = lineRespaced[i] - original_cam;
            ofVec3f normal = lineRespaced[i-1] - original_cam;
            
            normal.cross(v1);
            
            //if(normal.x <0)
                //normal = -normal;
            //else if(normal.x ==0 && normal.y<0)
                //normal = -normal;
            //else if(normal.x ==0 && normal.y==0 &&normal.z<0)
                //normal = -normal;
            normal.normalize();
            float d = normal.dot(original_cam);
            
            for(int j = 0; j<Tlist.size(); ++j)
            {
                ofxRayTriangleIntersection  rtIntersect;
                vector<FaceTri>             tris;
                vector<Ray>                 rays;
                
                
                FaceTri tri;
                tri.v0 = Tlist[j].p[0];
                tri.v1 = Tlist[j].p[1];
                tri.v2 = Tlist[j].p[2];
                tris.push_back(tri);
                
                Ray ray;
                ray.rayOrig = original_cam;
                ray.rayEnd = original_cam+2*(lineRespaced[i-1]-original_cam);
                rays.push_back(ray);
                ray.rayEnd = original_cam+2*(lineRespaced[i]-original_cam);
                rays.push_back(ray);
                
                rtIntersect.checkMeshIntersection(rays, tris);
                if(rtIntersect.intersecctInfos.size() == 2)
                {
                    ofPoint line_p[2];
                    line_p[0] = rtIntersect.intersecctInfos[0].intersectPos;
                    line_p[1] = rtIntersect.intersecctInfos[1].intersectPos;
                    Tlist[j].line_seg.push_back(line_p);
                    //break;
                    
                    vector<ofPoint> not_cut;
                    
                    //of_triangle test = Tlist[j];
                    if(normal.dot(Tlist[j].p[0]) > d)
                        add_cut_point(cut_point, Tlist[j].p[0]);
                    else
                        not_cut.push_back(Tlist[j].p[0]);
                    if(normal.dot(Tlist[j].p[1]) > d)
                        add_cut_point(cut_point, Tlist[j].p[1]);
                    else
                        not_cut.push_back(Tlist[j].p[1]);
                    if(normal.dot(Tlist[j].p[2]) > d)
                        add_cut_point(cut_point, Tlist[j].p[2]);
                    else
                        not_cut.push_back(Tlist[j].p[2]);
                    
                    
                    cut_line.push_back(line_p);
                    if(not_cut.size()==2)
                    {
                        ofPoint added_triangle[3];
                        added_triangle[0] = line_p[0];
                        added_triangle[1] = line_p[1];
                        added_triangle[2] = not_cut[0];
                        more_to_add.push_back(added_triangle);
                        more_to_add.back().normal[0] = Tlist[j].normal[0];
                        more_to_add.back().normal[1] = Tlist[j].normal[1];
                        more_to_add.back().normal[2] = Tlist[j].normal[2];
                        
                        added_triangle[0] = not_cut[0];
                        added_triangle[1] = not_cut[1];
                        added_triangle[2] = line_p[1];
                        more_to_add.push_back(added_triangle);
                        more_to_add.back().normal[0] = Tlist[j].normal[0];
                        more_to_add.back().normal[1] = Tlist[j].normal[1];
                        more_to_add.back().normal[2] = Tlist[j].normal[2];
                    }
                    else if(not_cut.size()==1)
                    {
                        ofPoint added_triangle[3];
                        added_triangle[0] = line_p[0];
                        added_triangle[1] = line_p[1];
                        added_triangle[2] = not_cut[0];
                        more_to_add.push_back(added_triangle);
                        more_to_add.back().normal[0] = Tlist[j].normal[0];
                        more_to_add.back().normal[1] = Tlist[j].normal[1];
                        more_to_add.back().normal[2] = Tlist[j].normal[2];
                    }
                    
                    which_intersect[j] = 0;
                    
                }
                else //if(rtIntersect.intersecctInfos.size() == 1)
                {
                    ofxRayTriangleIntersection  rtIntersect_1;
                    vector<FaceTri>             tris_1;
                    vector<Ray>                 rays_1;
                    
                    FaceTri tri;
                    tri.v0 = original_cam;
                    tri.v1 = original_cam+2*(lineRespaced[i-1]-original_cam);
                    tri.v2 = original_cam+2*(lineRespaced[i]-original_cam);
                    tris_1.push_back(tri);
                    
                    Ray ray;
                    ray.rayOrig = Tlist[j].p[0];
                    ray.rayEnd = Tlist[j].p[1];
                    rays_1.push_back(ray);
                    ray.rayEnd = Tlist[j].p[2];
                    rays_1.push_back(ray);
                    ray.rayOrig = Tlist[j].p[1];
                    rays_1.push_back(ray);
                    
                    rtIntersect_1.checkMeshIntersection(rays_1, tris_1);
                    if(rtIntersect_1.intersecctInfos.size() ==1 && rtIntersect.intersecctInfos.size() ==1)
                    {
                        ofPoint line_p[2];
                        line_p[0] = rtIntersect.intersecctInfos[0].intersectPos;
                        line_p[1] = rtIntersect_1.intersecctInfos[0].intersectPos;
                        Tlist[j].line_seg.push_back(line_p);
                        //break;
                        
                        vector<ofPoint> not_cut;
                        
                        //of_triangle test = Tlist[j];
                        
                        if(normal.dot(Tlist[j].p[0]) > d)
                            add_cut_point(cut_point, Tlist[j].p[0]);
                        else
                            not_cut.push_back(Tlist[j].p[0]);
                        if(normal.dot(Tlist[j].p[1]) > d)
                            add_cut_point(cut_point, Tlist[j].p[1]);
                        else
                            not_cut.push_back(Tlist[j].p[1]);
                        if(normal.dot(Tlist[j].p[2]) > d)
                            add_cut_point(cut_point, Tlist[j].p[2]);
                        else
                            not_cut.push_back(Tlist[j].p[2]);
                        
                        
                        cut_line.push_back(line_p);
                        if(not_cut.size()==2)
                        {
                            ofPoint added_triangle[3];
                            added_triangle[0] = line_p[0];
                            added_triangle[1] = line_p[1];
                            added_triangle[2] = not_cut[0];
                            more_to_add.push_back(added_triangle);
                            more_to_add.back().normal[0] = Tlist[j].normal[0];
                            more_to_add.back().normal[1] = Tlist[j].normal[1];
                            more_to_add.back().normal[2] = Tlist[j].normal[2];
                            
                            added_triangle[0] = not_cut[0];
                            added_triangle[1] = not_cut[1];
                            added_triangle[2] = line_p[1];
                            more_to_add.push_back(added_triangle);
                            more_to_add.back().normal[0] = Tlist[j].normal[0];
                            more_to_add.back().normal[1] = Tlist[j].normal[1];
                            more_to_add.back().normal[2] = Tlist[j].normal[2];
                        }
                        else if(not_cut.size()==1)
                        {
                            ofPoint added_triangle[3];
                            added_triangle[0] = line_p[0];
                            added_triangle[1] = line_p[1];
                            added_triangle[2] = not_cut[0];
                            more_to_add.push_back(added_triangle);
                            more_to_add.back().normal[0] = Tlist[j].normal[0];
                            more_to_add.back().normal[1] = Tlist[j].normal[1];
                            more_to_add.back().normal[2] = Tlist[j].normal[2];
                        }

                         
                        which_intersect[j] = 0;
                    }
                    else if(rtIntersect_1.intersecctInfos.size() ==2 )
                    {
                        ofPoint line_p[2];
                        line_p[0] = rtIntersect_1.intersecctInfos[0].intersectPos;
                        line_p[1] = rtIntersect_1.intersecctInfos[1].intersectPos;
                        Tlist[j].line_seg.push_back(line_p);
                        //break;
                        
                        vector<ofPoint> not_cut;
                        
                        //of_triangle test = Tlist[j];
                        if(normal.dot(Tlist[j].p[0]) > d)
                            add_cut_point(cut_point, Tlist[j].p[0]);
                        else
                            not_cut.push_back(Tlist[j].p[0]);
                        if(normal.dot(Tlist[j].p[1]) > d)
                            add_cut_point(cut_point, Tlist[j].p[1]);
                        else
                            not_cut.push_back(Tlist[j].p[1]);
                        if(normal.dot(Tlist[j].p[2]) > d)
                            add_cut_point(cut_point, Tlist[j].p[2]);
                        else
                            not_cut.push_back(Tlist[j].p[2]);
                        
                        cut_line.push_back(line_p);
                        if(not_cut.size()==2)
                        {
                            ofPoint added_triangle[3];
                            added_triangle[0] = line_p[0];
                            added_triangle[1] = line_p[1];
                            added_triangle[2] = not_cut[0];
                            more_to_add.push_back(added_triangle);
                            more_to_add.back().normal[0] = Tlist[j].normal[0];
                            more_to_add.back().normal[1] = Tlist[j].normal[1];
                            more_to_add.back().normal[2] = Tlist[j].normal[2];
                            
                            
                            added_triangle[0] = not_cut[0];
                            added_triangle[1] = not_cut[1];
                            added_triangle[2] = line_p[1];
                            more_to_add.push_back(added_triangle);
                            more_to_add.back().normal[0] = Tlist[j].normal[0];
                            more_to_add.back().normal[1] = Tlist[j].normal[1];
                            more_to_add.back().normal[2] = Tlist[j].normal[2];
                        }
                        else if(not_cut.size()==1)
                        {
                            ofPoint added_triangle[3];
                            added_triangle[0] = line_p[0];
                            added_triangle[1] = line_p[1];
                            added_triangle[2] = not_cut[0];
                            more_to_add.push_back(added_triangle);
                            more_to_add.back().normal[0] = Tlist[j].normal[0];
                            more_to_add.back().normal[1] = Tlist[j].normal[1];
                            more_to_add.back().normal[2] = Tlist[j].normal[2];
                        }

                         
                        which_intersect[j] = 0;
                    }
                }
            }
        //}
    
    vector<of_triangle> tmp_triangle;
    for (int i =0; i<Tlist.size() ; ++i)
    {
        if(which_intersect[i])
            tmp_triangle.push_back(Tlist[i]);
    }
    
    for (int i =0; i<more_to_add.size() ; ++i)
    {
        
        tmp_triangle.push_back(more_to_add[i]);
    }
    
    
    Tlist = tmp_triangle;
    T_num = Tlist.size();
    
    cut_construct(normal);
    
}

void testApp::paint_line()
{
    ofPolyline lineRespaced = line;
    //lineRespaced = lineRespaced.getResampledBySpacing(20);
    ofPoint original_cam = ofPoint(512,384,400);
    vector<of_triangle> tmp;
    vector<int> which_t;
    for(int j = 0; j<Tlist.size(); ++j)
        {
            if(Tlist[j].p[0].z>0 && Tlist[j].p[1].z>0 &&Tlist[j].p[2].z>0)
            {
                tmp.push_back(Tlist[j]);
                which_t.push_back(j);
            }
            //which_t.push_back(j);
        }
    //tmp = Tlist;
    for(int i =1; i <lineRespaced.size(); ++i)
    {
        ofVec3f v1 = lineRespaced[i-1] - original_cam;
        ofVec3f v2 = lineRespaced[i] - original_cam;
        v2.cross(v1);
        if(v2.y<0)
            v2 = -v2;
        v2.normalize();
        ofVec3f normal = lineRespaced[i] - lineRespaced[i-1];
        normal.cross(v2);
        
//        vector<int> tmp;
        if(normal.z < 0)
            normal = -normal;
        normal.normalize();
        float d = normal.dot(lineRespaced[i]);
        
        for(int j = 0; j<tmp.size(); ++j)
        {
            float angle = normal.angle(tmp[j].normal[0]);
            if( 1)//angle<=90)
            {
                ofxRayTriangleIntersection  rtIntersect;
                vector<FaceTri>             tris;
                vector<Ray>                 rays;
                
                
                FaceTri tri;
                tri.v0 = tmp[j].p[0];
                tri.v1 = tmp[j].p[1];
                tri.v2 = tmp[j].p[2];
                tris.push_back(tri);
                
                Ray ray;
                ray.rayOrig = original_cam;
                ray.rayEnd = lineRespaced[i-1];
                rays.push_back(ray);
                ray.rayEnd = lineRespaced[i];
                rays.push_back(ray);
                
                rtIntersect.checkMeshIntersection(rays, tris);
                if(rtIntersect.intersecctInfos.size() == 2 && rtIntersect.intersecctInfos[0].distance<1 && rtIntersect.intersecctInfos[1].distance<1)
                {
                    ofPoint line_p[2];
                    line_p[0] = rtIntersect.intersecctInfos[0].intersectPos;
                    line_p[1] = rtIntersect.intersecctInfos[1].intersectPos;
                    Tlist[which_t[j]].line_seg.push_back(line_p);
                    triInRing.push_back(which_t[j]);
                    //break;
                }
                else //if(rtIntersect.intersecctInfos.size() == 1)
                {
                    ofxRayTriangleIntersection  rtIntersect_1;
                    vector<FaceTri>             tris_1;
                    vector<Ray>                 rays_1;
                    
                    FaceTri tri;
                    tri.v0 = original_cam;
                    tri.v1 = lineRespaced[i-1];
                    tri.v2 = lineRespaced[i];
                    tris_1.push_back(tri);
                    
                    Ray ray;
                    ray.rayOrig = tmp[j].p[0];
                    ray.rayEnd = tmp[j].p[1];
                    rays_1.push_back(ray);
                    ray.rayEnd = tmp[j].p[2];
                    rays_1.push_back(ray);
                    ray.rayOrig = tmp[j].p[1];
                    rays_1.push_back(ray);

                    rtIntersect_1.checkMeshIntersection(rays_1, tris_1);
                    if(  rtIntersect_1.intersecctInfos.size() ==1 )//&& rtIntersect_1.intersecctInfos[0].distance<1 && rtIntersect.intersecctInfos[0].distance<1)
                    {
                        ofPoint line_p[2];
                        line_p[0] = rtIntersect.intersecctInfos[0].intersectPos;
                        line_p[1] = rtIntersect_1.intersecctInfos[0].intersectPos;
                        Tlist[which_t[j]].line_seg.push_back(line_p);
                        triInRing.push_back(which_t[j]);
                        //break;
                    }
                    /*
                }
                else
                {
                    ofxRayTriangleIntersection  rtIntersect_1;
                    vector<FaceTri>             tris_1;
                    vector<Ray>                 rays_1;
                    
                    FaceTri tri;
                    tri.v0 = original_cam;
                    tri.v1 = lineRespaced[i-1];
                    tri.v2 = lineRespaced[i];
                    tris_1.push_back(tri);
                    
                    Ray ray;
                    ray.rayOrig = tmp[j].p[0];
                    ray.rayEnd = tmp[j].p[1];
                    rays_1.push_back(ray);
                    ray.rayEnd = tmp[j].p[2];
                    rays_1.push_back(ray);
                    ray.rayOrig = tmp[j].p[1];
                    rays_1.push_back(ray);
                    
                    rtIntersect_1.checkMeshIntersection(rays_1, tris_1);*/
                    else if(rtIntersect_1.intersecctInfos.size() ==2 )//&& rtIntersect_1.intersecctInfos[0].distance<1 && rtIntersect_1.intersecctInfos[1].distance<1)
                    {
                        ofPoint line_p[2];
                        line_p[0] = rtIntersect_1.intersecctInfos[0].intersectPos;
                        line_p[1] = rtIntersect_1.intersecctInfos[1].intersectPos;
                        Tlist[which_t[j]].line_seg.push_back(line_p);
                        triInRing.push_back(which_t[j]);
                        //break;
                    }
                }
            }
        }
    }
    tmp.clear();
    line.clear();
}

void testApp::quarter_oval()
{
    vector<of_triangle> tmp_T;
    int tmp_T_num=0;
    for(int i =0; i<T_num; ++i)
    {
        ofPolyline poly1,poly2;
        poly1.clear();
        poly2.clear();
        vector<ofPoint> r_oval1,r_oval2;
        r_oval1.clear();
        r_oval2.clear();
        ofVec3f axis_x = ofVec3f(1,0,0);
        ofVec3f axis_z = ofVec3f(0,0,1);
        if(Tlist[i].counter[1])
        
        {
            
            ofPoint xy_point;
            xy_point = Tlist[i].p[0];
            xy_point.z = 0;
            ofVec3f v1 = Tlist[i].p[2] - xy_point;
            float rotate_angle = v1.angle(axis_x);
            ofPoint r_p1 = v1;
            r_p1.rotate(rotate_angle, axis_z);
            if(r_p1.angle(axis_x) > 1 )
            {
                rotate_angle = -rotate_angle;
                r_p1 = v1.rotate(rotate_angle, axis_z);
            }
            poly1.arc(xy_point,  r_p1.length(), Tlist[i].p[0].z,  270, 0, true, 20);
            for (int j =0; j<poly1.size(); ++j)
            {
                ofPoint tmp_collect = poly1[j];
                tmp_collect = tmp_collect.rotate(-90, xy_point, axis_x);
                tmp_collect = tmp_collect.rotate(-rotate_angle, xy_point, axis_z);
                r_oval1.push_back(tmp_collect);
            }
            
            xy_point = Tlist[i].p[1];
            xy_point.z = 0;
            v1 = Tlist[i].p[2] - xy_point;
            rotate_angle = v1.angle(axis_x);
            r_p1 = v1;
            r_p1.rotate(rotate_angle, axis_z);
            if(r_p1.angle(axis_x) > 1 )
            {
                rotate_angle = -rotate_angle;
                r_p1 = v1.rotate(rotate_angle, axis_z);
            }
            poly2.arc(xy_point, r_p1.length(), Tlist[i].p[1].z,  270, 0, true, 20);
            for (int j =0; j<poly2.size(); ++j)
            {
                ofPoint tmp_collect = poly2[j];
                tmp_collect = tmp_collect.rotate(-90, xy_point, axis_x);
                tmp_collect = tmp_collect.rotate(-rotate_angle, xy_point, axis_z);
                r_oval2.push_back(tmp_collect);
            }
        }
        else
        {
            ofPoint xy_point;
            xy_point = Tlist[i].p[0];
            xy_point.z = 0;
            ofVec3f v1 = Tlist[i].p[1] - xy_point;
            float rotate_angle = v1.angle(axis_x);
            ofPoint r_p1 = v1;
            r_p1.rotate(rotate_angle, axis_z);
            if(r_p1.angle(axis_x) > 1 )
            {
                rotate_angle = -rotate_angle;
                r_p1 = v1.rotate(rotate_angle, axis_z);
            }
            poly1.arc(xy_point,  r_p1.length(), Tlist[i].p[0].z,  270, 0, true,20);
            for (int j =0; j<poly1.size(); ++j)
            {
                ofPoint tmp_collect = poly1[j];
                tmp_collect = tmp_collect.rotate(-90, xy_point, axis_x);
                tmp_collect = tmp_collect.rotate(-rotate_angle, xy_point, axis_z);
                r_oval1.push_back(tmp_collect);
            }
            
            v1 = Tlist[i].p[2] - xy_point;
            rotate_angle = v1.angle(axis_x);
            r_p1 = v1;
            r_p1.rotate(rotate_angle, axis_z);
            if(r_p1.angle(axis_x) > 1 )
            {
                rotate_angle = -rotate_angle;
                r_p1 = v1.rotate(rotate_angle, axis_z);
            }
            poly2.arc(xy_point, r_p1.length(), Tlist[i].p[0].z,  270, 0, true, 20);
            for (int j =0; j<poly2.size(); ++j)
            {
                ofPoint tmp_collect = poly2[j];
                tmp_collect = tmp_collect.rotate(-90, xy_point, axis_x);
                tmp_collect = tmp_collect.rotate(-rotate_angle, xy_point, axis_z);
                r_oval2.push_back(tmp_collect);
            }
        }
        
        bool turn_1 = true;
        ofPoint T_p1,T_p2,T_p3;
        T_p1 = r_oval1[0];
        T_p3 = r_oval2[0];
        for(int j=1,k=1;j<r_oval1.size() || k<r_oval2.size();)
        {
            ofPoint T_tmppoint[3];
            if(turn_1)
            {
                T_p2 = T_p3;
                T_p3 = r_oval1[j++];
            }
            else
            {
                T_p1 = T_p3;
                T_p3 = r_oval2[k++];
            }
            
            turn_1 = !turn_1;
            T_tmppoint[0] = T_p1;
            T_tmppoint[1] = T_p2;
            T_tmppoint[2] = T_p3;
            tmp_T.push_back(T_tmppoint);
            
            ofVec3f v1 = T_p3 - T_p1;
            ofVec3f v2 = T_p2 - T_p1;
            v1 = v2.cross(v1);
            if(v1.z <0)
                v1 = -v1;
            v1 = v1.normalized();
            tmp_T[tmp_T_num].normal[0] = v1;
            tmp_T[tmp_T_num].normal[1] = v1;
            tmp_T[tmp_T_num].normal[2] = v1;
            ++tmp_T_num;
        }
        
    }
    Tlist.clear();
    Tlist = tmp_T;
    T_num = tmp_T_num;
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
                    Tlist[which_T[k]].p[which_P[k]].z = shift/1.5;
                    
                    ofVec3f tmp_n1 = Tlist[which_T[k]].p[2] - Tlist[which_T[k]].p[0];
                    ofVec3f tmp_n2 = Tlist[which_T[k]].p[1] - Tlist[which_T[k]].p[0];
                    tmp_n1 = tmp_n2.cross(tmp_n1);
                    if(tmp_n1.z < 0)
                        tmp_n1 = -tmp_n1;
                    //tmp_n1 = tmp_n1.normalized();
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
                            
                            //if(i == j)
                                //++i;
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
                            
                            //if(i == j)
                                //++i;
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
                            
                            //if(i == j)
                                //++i;
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

void testApp::clone()
{
    int T = T_num;
    for(int i = 0; i<T; ++i)
    {
        ofPoint tmp_p[3];
        tmp_p[0] = Tlist[i].p[0];
        tmp_p[1] = Tlist[i].p[1];
        tmp_p[2] = Tlist[i].p[2];
        of_triangle tmp = of_triangle(tmp_p);
        
        tmp.p[0].z = -tmp.p[0].z;
        tmp.p[1].z = -tmp.p[1].z;
        tmp.p[2].z = -tmp.p[2].z;
        
        ofVec3f v1 = tmp.p[1] -tmp.p[0];
        ofVec3f v2 = tmp.p[2] -tmp.p[0];
        v1 = v2.cross(v1);
        if(v1.z > 0)
            v1= -v1;
        v1.normalize();
        
        tmp.normal[0] = v1;
        tmp.normal[1] = v1;
        tmp.normal[2] = v1;
        
        ++T_num;
        Tlist.push_back(tmp);
    }
    
}


#pragma mark -
#pragma mark About OF

//--------------------------------------------------------------
void testApp::setup(){
    
    line.clear();
    Tlist.clear();
    T_num = 0;
    state = TO_CREATATION;
    
    ofSetSmoothLighting(true);
    pointLight.setDiffuseColor( ofFloatColor(1.f, 0.0f, 0.0f) );
    pointLight.setSpecularColor( ofFloatColor(1.0f, 0.0f, 0.0f));
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
    
    tmp_cam.resetTransform();
    tmp_cam.setFov(120);
    tmp_cam.clearParent();
    tmp_cam.setPosition(0, 0, 0);
    tmp_cam.lookAt(ofVec3f(512,384,0));
    
}

//--------------------------------------------------------------
void testApp::update(){
    //if(state == 4 && cut_point.size()!=0)
        //cut();
}

//--------------------------------------------------------------
void testApp::draw(){
    
    
    //ofBackgroundGradient(ofColor::white, ofColor(200,200,200), OF_GRADIENT_LINEAR);
    if(state == TO_CREATATION)
    {
        //tmp_cam.begin();
        cam.begin();
        line.draw();
        //cam.draw();
        //pointLight.draw();
        cam.end();
        //tmp_cam.end();
        
    }
    else if(state == CREATATAING)
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
    
    else if(state == TO_PAINT || state == TO_CUT || state == TO_REMOVE_TRI)
    {
        
        ofEnableDepthTest();
        
        ofEnableLighting();
        pointLight.enable();
        
        cam.begin();
        material.begin();
        ofSetColor(255,255,255);
        
        for(int i =0; i < T_num; ++i)
        {
            if (enableFace){
                glBegin(GL_TRIANGLES);
                glNormal3f(Tlist[i].normal[0].x,Tlist[i].normal[0].y,Tlist[i].normal[0].z);
                glVertex3f(Tlist[i].p[0].x,Tlist[i].p[0].y,Tlist[i].p[0].z);
                glNormal3f(Tlist[i].normal[1].x,Tlist[i].normal[1].y,Tlist[i].normal[1].z);
                glVertex3f(Tlist[i].p[1].x,Tlist[i].p[1].y,Tlist[i].p[1].z);
                glNormal3f(Tlist[i].normal[2].x,Tlist[i].normal[2].y,Tlist[i].normal[2].z);
                glVertex3f(Tlist[i].p[2].x,Tlist[i].p[2].y,Tlist[i].p[2].z);
                glEnd();
            }
            else{
                Tlist[i].draw_wireframe();
            }
        }
        
        material.end();
        cam.end();
        
        ofDisableLighting();
        cam.begin();
        
        /*
        ofSetColor(255, 255, 255);
        for(int i =0; i < T_num; ++i)
        {
            Tlist[i].draw_wireframe();
        }
        */

        
        ofSetColor(0, 0, 0);
        for(int i =0; i < T_num; ++i)
        {
            //Tlist[i].draw_wireframe();
            for (int j= 0; j<Tlist[i].line_seg.size(); ++j)
            {
                ofLine(Tlist[i].line_seg[j].p[0], Tlist[i].line_seg[j].p[1]);
            }
        }
        
        cam.end();
        ofDisableDepthTest();

    }
    
    else if(state == PAINTING || state == CUTTING || state == DRAWING_RING)
    {
        
        ofEnableDepthTest();
        
        ofEnableLighting();
        pointLight.enable();
        
        cam.begin();
        material.begin();
        ofSetColor(255,255,255);
        
        for(int i =0; i < T_num; ++i)
        {
            if (enableFace){
                glBegin(GL_TRIANGLES);
                glNormal3f(Tlist[i].normal[0].x,Tlist[i].normal[0].y,Tlist[i].normal[0].z);
                glVertex3f(Tlist[i].p[0].x,Tlist[i].p[0].y,Tlist[i].p[0].z);
                glNormal3f(Tlist[i].normal[1].x,Tlist[i].normal[1].y,Tlist[i].normal[1].z);
                glVertex3f(Tlist[i].p[1].x,Tlist[i].p[1].y,Tlist[i].p[1].z);
                glNormal3f(Tlist[i].normal[2].x,Tlist[i].normal[2].y,Tlist[i].normal[2].z);
                glVertex3f(Tlist[i].p[2].x,Tlist[i].p[2].y,Tlist[i].p[2].z);
                glEnd();
            }
            else{
                Tlist[i].draw_wireframe();
            }
        }
        
        material.end();
        cam.end();
        
        ofDisableLighting();
        
        
        cam.begin();
        
        ofSetColor(0, 0, 0);
        for(int i =0; i < T_num; ++i)
        {
            for (int j= 0; j<Tlist[i].line_seg.size(); ++j)
            {
                ofLine(Tlist[i].line_seg[j].p[0], Tlist[i].line_seg[j].p[1]);
            }
        }
        
        ofDisableDepthTest();
        
        ofSetColor(0, 0, 0);
        line.draw();
        
        cam.end();
    }

    else
    {
        ofEnableDepthTest();
        
        ofEnableLighting();
        pointLight.enable();
        
        cam.begin();
        material.begin();
        ofSetColor(255,255,255);
        
        for(int i =0; i < T_num; ++i)
        {
            if (enableFace){
                glBegin(GL_TRIANGLES);
                glNormal3f(Tlist[i].normal[0].x,Tlist[i].normal[0].y,Tlist[i].normal[0].z);
                glVertex3f(Tlist[i].p[0].x,Tlist[i].p[0].y,Tlist[i].p[0].z);
                glNormal3f(Tlist[i].normal[1].x,Tlist[i].normal[1].y,Tlist[i].normal[1].z);
                glVertex3f(Tlist[i].p[1].x,Tlist[i].p[1].y,Tlist[i].p[1].z);
                glNormal3f(Tlist[i].normal[2].x,Tlist[i].normal[2].y,Tlist[i].normal[2].z);
                glVertex3f(Tlist[i].p[2].x,Tlist[i].p[2].y,Tlist[i].p[2].z);
                glEnd();
            }
            else{
                Tlist[i].draw_wireframe();
            }
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

}

//--------------------------------------------------------------
void testApp::keyPressed(int key){
    //printf("%d", key);
    switch (key) {
        case 'h':
            draw_c = !draw_c;
            break;
            
        case 'c':
            state = TO_CUT;
            break;
            
        case 'p':
            state = TO_PAINT;
            break;
            
        case 's':
            state = TO_DRAW_RING;
            break;
        
        case 'z':
            enableFace = !enableFace;
            break;
            
        case 'f':
            if (state == TO_REMOVE_TRI) {
                state = TO_EXTRUDE;
                
                // create ring triangle
                // remove triangle in the ring
                rotate(-90, ofVec3f(1,0,0));
                
                // construct ring shape
                createRing();
                //removeTriInRing();
                // remove origin triangle the ring pass
                removeRing();
            }
            break;
            
        case 'o':
            for(int i = 0; i<Tlist.size(); ++i)
                Tlist[i].line_seg.clear();
            break;
            
        case 'e':
            prune_1();
            prune_1(); //has bug : why must do prune_1 twice?
            prune_2();
            
            elevate();
            quarter_oval();
            cast();
            clone();
            //elevated_T = true;
            state = TO_PAINT;
            line.clear();
            break;
        
        case 'r':
            line.clear();
            Tlist.clear();
            T_num = 0;
            state = TO_CREATATION;
            break;
            
        case 'q':
            cut();
            break;
            
        case OF_KEY_UP:
            translate(10, ofVec3f(0, 0, 1));
            break;
                      
        case OF_KEY_DOWN:
            translate(10, ofVec3f(0, 0, -1));
            break;
            
        default:
            break;
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
    if (button == 2) {
        rotate( (y - old_y) * 0.3 , ofVec3f(1, 0, 0));
        rotate( (x - old_x) * 0.3 , ofVec3f(0, 1, 0));
        old_x = x;
        old_y = y;
        return;
    }
    else if (button == 1){
        translate((x - old_x) * 1, ofVec3f(1, 0, 0));
        translate((y - old_y) * 1, ofVec3f(0, -1, 0));
        old_x = x;
        old_y = y;
        return;
    }
    
    if(state == PAINTING)
        line.addVertex(ofPoint(x, 768-y));
    else if (state == CUTTING)
        line.addVertex(ofPoint(x, 768-y));
    else
        line.addVertex(ofPoint(x, 768-y));
}

//--------------------------------------------------------------
void testApp::mousePressed(int x, int y, int button){
    if (button == 0) { // left
        if(state == TO_CREATATION)
            line.addVertex(ofPoint(x, 768-y));
        else if(state == TO_PAINT)
        {
            line.clear();
            line.addVertex(ofPoint(x, 768-y));
            state = PAINTING;
        }
        else if(state == TO_CUT)
        {
            line.clear();
            line.addVertex(ofPoint(x, 768-y));
            state = CUTTING;
        }
        else if(state == TO_DRAW_RING){
            line.clear();
            line.addVertex(ofPoint(x, 768-y));
            state = DRAWING_RING;
            //printf("(%d, %d)\n", x, y);
        }
    }
    else if (button == 1){ // middle
        stateSave = state;
        old_x = x;
        old_y = y;
        state = TRANSLATE;
    }
    else if (button == 2){ // right
        stateSave = state;
        old_x = x;
        old_y = y;
        state = ROTATE;
    }
}

//--------------------------------------------------------------
void testApp::mouseReleased(int x, int y, int button){
    if (button == 2) {
        state = stateSave;
        return;
    }
    if (button == 1){
        state = stateSave;
        return;
    }
    
    if(state == TO_CREATATION)
    {
        if (line.size() > 2){
        
        release = true;
        state = CREATATAING;
        
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
    else if(state == PAINTING)
    {
        state = TO_PAINT;
        paint_line();
    }
    else if(state == CUTTING)
    {
        state = TO_CUT;
        cut_plane();
        cut();
    }
    else if(state == DRAWING_RING){
        state = TO_REMOVE_TRI;
        
        // ring in surface
        line.addVertex(line[0]);
        paint_line();
        

        
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

#pragma mark -
#pragma mark About extrusion
void testApp::createRing(){
    
    // reconstruct the triangle along the ring
    for (int i = 0; i < triInRing.size(); i++) {
        of_triangle tri = Tlist[triInRing[i]];
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

void testApp::removeTriInRing(){
    
    float level = Tlist[triInRing[0]].line_seg[0].p[0].y;
    
    // remove the triangle in the ring
    // There would be some bug is the shape is not convex
    for (vector<of_triangle>::iterator it = Tlist.begin(); it != Tlist.end();) {
        if ((*it).p[0].y > level && (*it).p[1].y > level && (*it).p[2].y > level)
            it = Tlist.erase(it);
        else
            it++;
    }
    
    triInRing.clear();
}

void testApp::removeRing(){
    std::sort(triInRing.begin(), triInRing.end());
    triInRing.erase(unique(triInRing.begin(), triInRing.end()), triInRing.end());
    for (int i = triInRing.size()-1; i >= 0; i--) {
        Tlist.erase(Tlist.begin() + triInRing[i]);
    }
}


#pragma mark operating 
void testApp::rotate(float theta, ofVec3f dir){
    for(int i = 0; i<T_num; ++i)
    {
        Tlist[i].p[0] = Tlist[i].p[0].rotate(theta, ofVec3f(512,384,0), dir);
        Tlist[i].p[1] = Tlist[i].p[1].rotate(theta, ofVec3f(512,384,0), dir);
        Tlist[i].p[2] = Tlist[i].p[2].rotate(theta, ofVec3f(512,384,0), dir);
        
        Tlist[i].normal[0] = Tlist[i].normal[0].rotate(theta, dir);
        Tlist[i].normal[1] = Tlist[i].normal[1].rotate(theta, dir);
        Tlist[i].normal[2] = Tlist[i].normal[2].rotate(theta, dir);
        
        for(int j = 0; j < Tlist[i].line_seg.size(); ++j)
        {
            Tlist[i].line_seg[j].p[0]=Tlist[i].line_seg[j].p[0].rotate(theta, ofVec3f(512,384,0), dir);
            Tlist[i].line_seg[j].p[1]=Tlist[i].line_seg[j].p[1].rotate(theta, ofVec3f(512,384,0), dir);
        }
        
    }
}



void testApp::translate(float dist, ofVec3f dir){
    for(int i = 0; i<T_num; ++i)
    {
        Tlist[i].p[0] = Tlist[i].p[0] + dist * dir;
        Tlist[i].p[1] = Tlist[i].p[1] + dist * dir;
        Tlist[i].p[2] = Tlist[i].p[2] + dist * dir;
        
        for(int j = 0; j < Tlist[i].line_seg.size(); ++j)
        {
            Tlist[i].line_seg[j].p[0]=Tlist[i].line_seg[j].p[0] + dir;
            Tlist[i].line_seg[j].p[1]=Tlist[i].line_seg[j].p[1] + dir;
        }
        
    }
}