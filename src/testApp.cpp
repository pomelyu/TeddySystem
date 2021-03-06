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
                    triBelongToRing.push_back(which_t[j]);
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
                        triBelongToRing.push_back(which_t[j]);
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
                        triBelongToRing.push_back(which_t[j]);
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
        
        /*
        // ========= for debug ================ //
        if (isTest){
            for (int i = 0; i < tmp1.size(); i++) {
                ofDrawSphere(Tlist[tmp1[i]].p[tmp2[i]], 2);
            }
            ofSetColor(0, 0, 0);
            for (int i = 0; i < triBelongToRing.size(); i++) {
                ofDrawSphere(Tlist[triBelongToRing[i]].line_seg.front().p[0], 0.1);
                ofDrawSphere(Tlist[triBelongToRing[i]].line_seg.back().p[1], 0.1);
            }
        }
        // ===================================== //
        */

        for(int i =0; i < T_num; ++i)
        {
            if (enableFace){
                
                // ========= for debug ================ //
                if (isTest) {
                    if (Tlist[i].tColor == NONE)
                        ofSetColor(0,0,255);
                    else if (Tlist[i].tColor == GRAY)
                        ofSetColor(0, 255, 0);
                    else if (Tlist[i].tColor == BLACK)
                        ofSetColor(255, 0, 0);
                    else
                        ofSetColor(255, 255, 255);
                }
                // ===================================== //
                
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
        
        // For test
        ofSetColor(0, 0, 0);
        for (int i = 0; i < triBelongToRing.size(); i++) {
            of_triangle* tri = &(Tlist[triBelongToRing[i]]);
            for (int j = 0 ; j < tri->line_seg.size(); j++) {
                ofLine(tri->line_seg[j].p[0], tri->line_seg[j].p[1]);
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
                //rotate(-90, ofVec3f(1,0,0));
                // sorting the line edge in each triangle
                sortLineSeg();
                simplifyLine();
                // construct ring shape
                seperateTri();
                printf( "%lu\n", triInsideRing.size());
            }
            break;
            
        case 'x':
            if (state == TO_EXTRUDE) {
                createRing();
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
            triBelongToRing.clear();
            triInsideRing.clear();
            colorInsideRing = NONE;
            isTest = false;
            break;
            
        case 'q':
            cut();
            break;
        
        case 't':
            isTest = !isTest;
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
            Tlist[i].line_seg[j].p[0]=Tlist[i].line_seg[j].p[0] + dist * dir;
            Tlist[i].line_seg[j].p[1]=Tlist[i].line_seg[j].p[1] + dist * dir;
        }
        
    }
}

#pragma mark -
#pragma mark About extrusion

void testApp::sortLineSeg(){
    
    // remove duplicate in triBelongToRing
    std::sort(triBelongToRing.begin(), triBelongToRing.end());
    triBelongToRing.erase(unique(triBelongToRing.begin(), triBelongToRing.end()), triBelongToRing.end());
    
    for (int i = 0; i < triBelongToRing.size(); i++) {
        vector<of_edge> edge = Tlist[triBelongToRing[i]].line_seg;
        vector<of_edge> sortedLine;
        
        int n = 0;
        while (edge.size() > 0) {
            vector<of_edge> vec;
            
            vec.push_back(edge[0]);
            edge.erase(edge.begin());
            
            ofPoint head;
            ofPoint tail;
            
            while (edge.size() > 0) {
                head = vec.front().p[0];
                tail = vec.back().p[1];
                
                bool isFound = false;
                for (int j = 0; j < edge.size(); j++) {
                    if (edge[j].p[0] == head) {
                        ofPoint tmp = edge[j].p[0];
                        edge[j].p[0] = edge[j].p[1];
                        edge[j].p[1] = tmp;
                        vec.insert(vec.begin(), edge[j]);
                        edge.erase(edge.begin() + j);
                        isFound = true;
                        break;
                    }
                    else if (edge[j].p[1] == head){
                        vec.insert(vec.begin(), edge[j]);
                        edge.erase(edge.begin() + j);
                        isFound = true;
                        break;
                    }
                    else if (edge[j].p[0] == tail){
                        vec.push_back(edge[j]);
                        edge.erase(edge.begin() + j);
                        isFound = true;
                        break;
                    }
                    else if (edge[j].p[1] == tail){
                        ofPoint tmp = edge[j].p[0];
                        edge[j].p[0] = edge[j].p[1];
                        edge[j].p[1] = tmp;
                        vec.push_back(edge[j]);
                        edge.erase(edge.begin() + j);
                        isFound = true;
                        break;
                    }
                }
                
                if (!isFound)
                    break;
            }
            for (int j = 0 ; j < vec.size(); j++) {
                sortedLine.push_back(vec[j]);
            }
            if (++n == 1) {
                Tlist[triBelongToRing[i]].numSplit = vec.size();
            }
        }
        Tlist[triBelongToRing[i]].line_seg = sortedLine;
        
        if ( n > 1)
            printf("line with %d segment\n", n);
    }
}


void testApp::seperateTri(){
    
    // Draw all triangle to NONE
    for (int i = 0; i < Tlist.size(); i++)
        Tlist[i].tColor = NONE;
    
    // Draw all triangle in ring to WHITE
    for (int i = 0; i < triBelongToRing.size(); i++)
        Tlist[triBelongToRing[i]].tColor =WHITE;
    
    // Draw one of vertex in a ring triangle to BLACK
    // Draw the triangle to GRAY
    // GRAY means the vertex of triangle has been drawed
    of_triangle* tri = &(Tlist[triBelongToRing[0]]);
    
    drawRingTriangle(tri, 0);
    
    queue<of_triangle*> triStack;
    triStack.push(tri);
    ofPoint pt;
    
    // By BFS
    while (triStack.size() != 0) {
        tri = triStack.front();
        triStack.pop();
        
        // get the none BLACK vertex
        for (int idx = 0; idx < 3; idx++) {
            if (tri->vColor[idx] == GRAY) {
                tri->vColor[idx] == BLACK;
                pt = tri->p[idx];
                
                // Who has this point
                for (int i = 0; i < Tlist.size(); i++) {
                    of_triangle* tr = &(Tlist[i]);
                    
                    for (int sameIdx = 0; sameIdx < 3; sameIdx++) {
                        if (check_point_same(pt, tr->p[sameIdx])){
                            if (tr->tColor == WHITE) {
                                drawRingTriangle(tr, sameIdx);
                                triStack.push(tr);
                            }
                            else if (tr->tColor == NONE){
                                for (int k = 0; k < 3; k++)
                                    tr->vColor[k] = GRAY;
                                tr->vColor[sameIdx] = BLACK;
                                tr->tColor = BLACK;
                                triStack.push(tr);
                            }
                            else
                                tr->vColor[sameIdx] = BLACK;
                            break;
                        }
                    }
                }
            }
        }
    }
    
    // Choose which color is in Ring
    vector<int> vec_BLACK;
    vector<int> vec_NONE;
    
    for (int i = 0; i < Tlist.size(); i++) {
        if (Tlist[i].tColor == BLACK)
            vec_BLACK.push_back(i);
        else if (Tlist[i].tColor == NONE)
            vec_NONE.push_back(i);
        
    }
    
    if (vec_BLACK.size() > vec_NONE.size()) {
        colorInsideRing = NONE;
        triInsideRing = vec_NONE;
    }
    else{
        colorInsideRing = BLACK;
        triInsideRing = vec_BLACK;
    }
}

void testApp::drawRingTriangle(of_triangle* tri, int baseIdx){
    
    //======== Just one line segment =======
    if (tri->numSplit == tri->line_seg.size()) {
        // Check if three vertex are in or out of ring
        bool threeVertexOutOfRing = false;
        
        for (int j = 0; j < 3; j++) {
            ofVec3f vec1 = tri->p[j] - tri->line_seg.front().p[0];
            ofVec3f vec2 = tri->p[j] - tri->line_seg.back().p[1];
            ofVec3f vec3 = tri->line_seg.front().p[0] - tri->line_seg.back().p[1];
            if (vec1.length() > 5 && vec2.length() > 5) {
                if (abs(1 - abs(vec1.dot(vec3) / ( vec1.length() * vec3.length())) ) < 0.001 &&
                    abs(1 - (vec1.dot(vec2)) / (vec1.length() * vec2.length())) < 0.001) {
                    for (int k = 0; k < 3; k++)
                        tri->vColor[k] = GRAY;
                    threeVertexOutOfRing = true;
                    
                    //tri->tColor = GRAY;
                    tri->tColor = YELLOW;
                    printf("All in ring. case 1\n");
                    break;
                }
            }
        }
        
        
        // Draw appropriate color to vertex
        if (!threeVertexOutOfRing) {
            ofVec3f norm = (tri->normal[0]).cross(tri->line_seg.front().p[0] - tri->line_seg.back().p[1]);
            float onePt = norm.dot(tri->p[baseIdx] - tri->line_seg.front().p[0]);
            for (int k = 0; k < 3; k++) {
                if (norm.dot(tri->p[k] - tri->line_seg.front().p[0]) * onePt > 0)
                    tri->vColor[k] = GRAY;
            }
            tri->tColor = GRAY;
        }
        
        tri->vColor[baseIdx] = BLACK;
    }
    //======== Two line segment =======
    else{
        int split = tri->numSplit;
        ofVec3f norm1 = (tri->normal[0]).cross(tri->line_seg.front().p[0] - tri->line_seg[split-1].p[1]);
        ofVec3f norm2 = (tri->normal[0]).cross(tri->line_seg[split].p[0] - tri->line_seg.back().p[1]);
        float base = norm1.dot(tri->p[baseIdx] - tri->line_seg.front().p[0]);
        base = base * norm2.dot(tri->p[baseIdx] - tri->line_seg.back().p[1]);
        
        int grayPt = 1;
        for (int k = 0; k < 3; k++) {
            float compare = norm1.dot(tri->p[k] - tri->line_seg.front().p[0]);
            compare = compare * norm2.dot(tri->p[k] - tri->line_seg.back().p[1]);
            
            if (compare * base > 0){
                tri->vColor[k] = GRAY;
                grayPt++;
            }
        }
        tri->vColor[baseIdx] = BLACK;
        if (grayPt == 3){
            tri->tColor = YELLOW;
            printf("All in ring. case 2\n");
        }
        else
            tri->tColor = GRAY;
    }
}


void testApp::simplifyLine(){
    for (int i = 0; i < triBelongToRing.size(); i++) {
        of_triangle* tri = &(Tlist[triBelongToRing[i]]);
        
        // ==== Only one line ====
        if (tri->numSplit == tri->line_seg.size()) {
            if (tri->line_seg.size() > 1){
                tri->line_seg = simplifyOneLine(tri, 0, tri->line_seg.size()-1);
                tri->numSplit = tri->line_seg.size();
            }
        }
        // ==== If two line ======
        else{
            vector<of_edge> tmp1 = simplifyOneLine(tri, 0, tri->numSplit-1);
            vector<of_edge> tmp2 = simplifyOneLine(tri, tri->numSplit, tri->line_seg.size()-1);
            tri->line_seg = tmp1;
            tri->numSplit = tmp1.size();
            for (int j = 0; j < tmp2.size(); j++) {
                tri->line_seg.push_back(tmp2[j]);
            }
        }
    }
}


vector<of_edge> testApp::simplifyOneLine(of_triangle* tri, int first, int last){
    
    ofPoint head = tri->line_seg[first].p[0];
    ofPoint tail = tri->line_seg[last].p[1];
    
    vector<ofPoint> vecLine;
    vecLine.push_back(head);
    
    for (int i = first + 1; i <= last; i++) {
        ofVec3f vec1 = tri->line_seg[i].p[0] - head;
        ofVec3f vec2 = tri->line_seg[i].p[0] - tail;
        
        float product = -1;
        if (vec1.length() > 1 && vec2.length() > 1)
            product = vec1.dot(vec2) / (vec1.length() * vec2.length());
        
        // if angle < 120 record it
        if (product > -0.5) {
            head = tri->line_seg[i].p[0];
            vecLine.push_back(head);
        }
    }
    vecLine.push_back(tail);
    vector<of_edge> vecEdge;
    for (int i = 0; i < vecLine.size()-1; i++) {
        vecEdge.push_back(of_edge(vecLine[i], vecLine[i+1]));
    }
    return vecEdge;
}

void testApp::createRing(){
    vector<int> triDelete;
    
    for (int i = 0; i < triBelongToRing.size(); i++) {
        triDelete.push_back(triBelongToRing[i]);
        of_triangle* tri = &(Tlist[triBelongToRing[i]]);
        
        // calculate how many vertex not in the ring
        vector<int> vOutsideRing;
        for (int j = 0; j < 3; j++) {
            if (tri->vColor[j] != colorInsideRing)
                vOutsideRing.push_back(j);
        }
        
        // ==== Just one line segment ====
        if (tri->numSplit == tri->line_seg.size()) {
            // If one vetex not in the ring, just connect the
            // vertex and each line segment to create new triangle
            if (vOutsideRing.size() == 1) {
                for (int j = 0; j < tri->line_seg.size(); j++) {
                    of_triangle tmpTri(tri->p[vOutsideRing[0]],
                                       tri->line_seg[j].p[0],
                                       tri->line_seg[j].p[1]);
                    tmpTri.copyNormal(tri->normal);
                    Tlist.push_back(tmpTri);
                    T_num++;
                }
            }
            // if two vertex not in the ring. calculate the lowest
            // point in line segment. Connect each line segment to
            // appropriate vertex to create new triangle
            else if (vOutsideRing.size() == 2){
                
                // Find the lowest point in line segment
                float value = 9999;
                int lowest = 0;
                for (int j = 0; j < tri->line_seg.size(); j++) {
                    ofVec3f v1 = tri->line_seg[j].p[0] - tri->p[vOutsideRing[0]];
                    ofVec3f v2 = tri->line_seg[j].p[0] - tri->p[vOutsideRing[1]];
                    float v = v1.dot(v2) / ( v1.length()*v2.length());
                    if (v < value) {
                        value = v;
                        lowest = j;
                    }
                }
                
                // connect the lowest line segment and two vertex
                of_triangle tmpTri(tri->line_seg[lowest].p[0],
                                   tri->p[vOutsideRing[0]],
                                   tri->p[vOutsideRing[1]]);
                tmpTri.copyNormal(tri->normal);
                Tlist.push_back(tmpTri);
                T_num++;
                
                // connect each line segment to appropriate vertex
                for (int j = 0; j < 3; j++) {
                    if (tri->vColor[j] == colorInsideRing) {
                        ofVec3f v1 = tri->p[j] - tri->p[vOutsideRing[0]];
                        ofVec3f v2 = tri->p[j] - tri->line_seg[0].p[0];
                        float parallel = v1.dot(v2) / ( v1.length()*v2.length() );
                        
                        // tri->p[j] -- line_seg[0].p[0] -- tri->p[vOutsideRing[0]] in the same line
                        if (abs(1 - abs(parallel)) < 0.001) {
                            
                            // Triangle before "lowest"
                            for (int k = 0 ; k < lowest; k++) {
                                of_triangle tmpTri(tri->p[vOutsideRing[0]],
                                                   tri->line_seg[k].p[0],
                                                   tri->line_seg[k].p[1]);
                                tmpTri.copyNormal(tri->normal);
                                Tlist.push_back(tmpTri);
                                T_num++;
                            }
                            
                            // Triangle after "lowest"
                            for (int k = lowest; k < tri->line_seg.size(); k++) {
                                of_triangle tmpTri(tri->p[vOutsideRing[1]],
                                                   tri->line_seg[k].p[0],
                                                   tri->line_seg[k].p[1]);
                                tmpTri.copyNormal(tri->normal);
                                Tlist.push_back(tmpTri);
                                T_num++;
                            }
                        }
                        // tri->p[j] -- line_seg[0].p[1] -- tri->p[vOutsideRing[0]] in the same line
                        else{
                            
                            // Triangle before "lowest"
                            for (int k = 0 ; k < lowest; k++) {
                                of_triangle tmpTri(tri->p[vOutsideRing[1]],
                                                   tri->line_seg[k].p[0],
                                                   tri->line_seg[k].p[1]);
                                tmpTri.copyNormal(tri->normal);
                                Tlist.push_back(tmpTri);
                                T_num++;
                            }
                            
                            // Triangle after "lowest"
                            for (int k = lowest; k < tri->line_seg.size(); k++) {
                                of_triangle tmpTri(tri->p[vOutsideRing[0]],
                                                   tri->line_seg[k].p[0],
                                                   tri->line_seg[k].p[1]);
                                tmpTri.copyNormal(tri->normal);
                                Tlist.push_back(tmpTri);
                                T_num++;
                            }
                        }
                    }
                }
                
            }
            // If three vertex out of ring
            // copy the origin triangle and do nothing
            else{
                of_triangle tmpTri(*tri);
                tmpTri.line_seg.clear();
                Tlist.push_back(tmpTri);
                T_num++;
            }
        }
        // ==== Two line segment ====
        else{
            // WTF
            if (vOutsideRing.size() == 1){
            }
            // Two vertex are outside of ring
            // connect each vertex to appropriate line segment
            else if (vOutsideRing.size() == 2){
                int idx = 0;
                for (; idx < 3; idx++) {
                    if (tri->vColor[idx] == colorInsideRing) {
                        break;
                    }
                }
                ofVec3f side1 = tri->p[vOutsideRing[0]] - tri->p[idx];
                ofVec3f side2 = tri->p[vOutsideRing[0]] - tri->p[vOutsideRing[1]];
                ofVec3f line1 = tri->p[vOutsideRing[0]] - tri->line_seg[0].p[0];
                ofVec3f line2 = tri->p[vOutsideRing[0]] - tri->line_seg[tri->numSplit-1].p[1];
                
                if (abs( side1.dot(side2)/(side1.length()*side2.length()) -
                         line1.dot(line2)/(line1.length()*line2.length()) ) < 0.0001) {
                    for (int j = 0; j < tri->numSplit; j++) {
                        of_triangle tmpTri(tri->p[vOutsideRing[0]],
                                           tri->line_seg[j].p[0],
                                           tri->line_seg[j].p[1]);
                        tmpTri.copyNormal(tri->normal);
                        Tlist.push_back(tmpTri);
                        T_num++;
                    }
                    for (int j = tri->numSplit; j < tri->line_seg.size(); j++) {
                        of_triangle tmpTri(tri->p[vOutsideRing[1]],
                                           tri->line_seg[j].p[0],
                                           tri->line_seg[j].p[1]);
                        tmpTri.copyNormal(tri->normal);
                        Tlist.push_back(tmpTri);
                        T_num++;
                    }
                }
                else{
                    for (int j = 0; j < tri->numSplit; j++) {
                        of_triangle tmpTri(tri->p[vOutsideRing[1]],
                                           tri->line_seg[j].p[0],
                                           tri->line_seg[j].p[1]);
                        tmpTri.copyNormal(tri->normal);
                        Tlist.push_back(tmpTri);
                        T_num++;
                    }
                    for (int j = tri->numSplit; j < tri->line_seg.size(); j++) {
                        of_triangle tmpTri(tri->p[vOutsideRing[0]],
                                           tri->line_seg[j].p[0],
                                           tri->line_seg[j].p[1]);
                        tmpTri.copyNormal(tri->normal);
                        Tlist.push_back(tmpTri);
                        T_num++;
                    }
                }
                
            }
            // Three vertex out of ring
            else{
            }
        }
        
    }
    
    for (int i = 0; i < triInsideRing.size(); i++)
        triDelete.push_back(triInsideRing[i]);
    
    std::sort(triDelete.begin(), triDelete.end());
    for (int i = triDelete.size()-1; i >= 0; i--) {
        Tlist.erase(Tlist.begin() + triDelete[i]);
    }
    
    triBelongToRing.clear();
    triInsideRing.clear();
}