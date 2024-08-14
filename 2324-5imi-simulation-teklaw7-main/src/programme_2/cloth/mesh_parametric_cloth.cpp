/*
**    TP CPE Lyon
**    Copyright (C) 2015 Damien Rohmer
**
**    This program is free software: you can redistribute it and/or modify
**    it under the terms of the GNU General Public License as published by
**    the Free Software Foundation, either version 3 of the License, or
**    (at your option) any later version.
**
**   This program is distributed in the hope that it will be useful,
**    but WITHOUT ANY WARRANTY; without even the implied warranty of
**    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
**    GNU General Public License for more details.
**
**    You should have received a copy of the GNU General Public License
**    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "mesh_parametric_cloth.hpp"

#include "../lib/common/error_handling.hpp"
#include <cmath>

namespace cpe
{


void mesh_parametric_cloth::update_force()
{

    int const Nu = size_u();
    int const Nv = size_v();
    // std::cout<<"Nu = "<<Nu<<std::endl;
    // std::cout<<"Nv = "<<Nv<<std::endl;
    int const N_total = Nu*Nv;
    ASSERT_CPE(static_cast<int>(force_data.size()) == Nu*Nv , "Error of size");

    float Kstruct = 10.0f;  //structural
    float Kshear = 8.0f;    //shearing
    float Kbend = 8.0f; //bending
    float L0 = 1/(Nu-1.0f); //length of a spring
    float Kw = 0.01f;   //wind force
    vec3 const wind_direction = vec3(1.0f,0.0f,0.0f);
    //Gravity
    static vec3 const g (0.0f,0.0f,-9.81f);
    vec3 const g_normalized = g/N_total;
    for(int ku=0 ; ku<Nu ; ++ku)
    {
        for(int kv=0 ; kv<Nv ; ++kv)
        {
            force(ku,kv) = g_normalized;    //gravity

            //spring structurels
            force(ku,kv) += spring_compute_force(ku,kv,ku+1,kv,L0,Kstruct);
            force(ku,kv) += spring_compute_force(ku,kv,ku,kv+1,L0,Kstruct);
            force(ku,kv) += spring_compute_force(ku,kv,ku-1,kv,L0,Kstruct);
            force(ku,kv) += spring_compute_force(ku,kv,ku,kv-1,L0,Kstruct);

            //spring shearing
            force(ku,kv) += spring_compute_force(ku,kv,ku+1,kv+1,L0*sqrt(2),Kshear);
            force(ku,kv) += spring_compute_force(ku,kv,ku-1,kv+1,L0*sqrt(2),Kshear);
            force(ku,kv) += spring_compute_force(ku,kv,ku+1,kv-1,L0*sqrt(2),Kshear);
            force(ku,kv) += spring_compute_force(ku,kv,ku-1,kv-1,L0*sqrt(2),Kshear);

            //spring bending
            force(ku,kv) += spring_compute_force(ku,kv,ku+2,kv,L0*2,Kbend);
            force(ku,kv) += spring_compute_force(ku,kv,ku,kv+2,L0*2,Kbend);
            force(ku,kv) += spring_compute_force(ku,kv,ku-2,kv,L0*2,Kbend);
            force(ku,kv) += spring_compute_force(ku,kv,ku,kv-2,L0*2,Kbend);

            //compute wind force
            force(ku,kv) += Kw*dot(wind_direction,normal(ku,kv))*normal(ku,kv); //wind force
        }
    }

    //*************************************************************//
    // TO DO, Calculer les forces s'appliquant sur chaque sommet
    //*************************************************************//
    //
    // ...
    //
    //
    //
    //*************************************************************//
    //Points fixes
     force(0,0) = vec3(0,0,0);
     force(0,Nv-1) = vec3(0,0,0);
    //  force(Nu-1,0) = vec3(0,0,0);
    //  force(Nu-1,Nv-1) = vec3(0,0,0);
}

void mesh_parametric_cloth::integration_step(float const dt)
{
    ASSERT_CPE(speed_data.size() == force_data.size(),"Incorrect size");
    ASSERT_CPE(static_cast<int>(speed_data.size()) == size_vertex(),"Incorrect size");


    int const Nu = size_u();
    int const Nv = size_v();
    float const mu = 0.1f;
    //*************************************************************//
    // TO DO: Calculer l'integration numerique des positions au cours de l'intervalle de temps dt.
    //*************************************************************//
    //
    // ...
    //
    //
    //
    //*************************************************************//

    for(int ku=0; ku<Nu ; ++ku)
    {
        for(int kv=0 ; kv<Nv ; ++kv)
        {
            vec3& p = vertex(ku,kv);
            vec3& v = speed(ku,kv);
            v = (1-mu*dt)*v + dt*force(ku,kv);
            p += dt*v;


            // intersection avec le plan
            intersect_plane(ku,kv);
            //intersection avec la sphère
            intersect_sphere(ku,kv,{0.5f,0.05f,-1.1f},0.2f);
        }
    }




    //security check (throw exception if divergence is detected)
    static float const LIMIT=30.0f;
    for(int ku=0 ; ku<Nu ; ++ku)
    {
        for(int kv=0 ; kv<Nv ; ++kv)
        {
            vec3 const& p = vertex(ku,kv);

            if( norm(p) > LIMIT )
            {
                throw exception_divergence("Divergence of the system",EXCEPTION_PARAMETERS_CPE);
            }
        }
    }

}

void mesh_parametric_cloth::set_plane_xy_unit(int const size_u_param,int const size_v_param)
{
    mesh_parametric::set_plane_xy_unit(size_u_param,size_v_param);

    int const N = size_u()*size_v();
    speed_data.resize(N);
    force_data.resize(N);
}

vec3 const& mesh_parametric_cloth::speed(int const ku,int const kv) const
{
    ASSERT_CPE(ku >= 0 , "Value ku ("+std::to_string(ku)+") should be >=0 ");
    ASSERT_CPE(ku < size_u() , "Value ku ("+std::to_string(ku)+") should be < size_u ("+std::to_string(size_u())+")");
    ASSERT_CPE(kv >= 0 , "Value kv ("+std::to_string(kv)+") should be >=0 ");
    ASSERT_CPE(kv < size_v() , "Value kv ("+std::to_string(kv)+") should be < size_v ("+std::to_string(size_v())+")");

    int const offset = ku + size_u()*kv;

    ASSERT_CPE(offset < static_cast<int>(speed_data.size()),"Internal error");

    return speed_data[offset];
}

vec3& mesh_parametric_cloth::speed(int const ku,int const kv)
{
    ASSERT_CPE(ku >= 0 , "Value ku ("+std::to_string(ku)+") should be >=0 ");
    ASSERT_CPE(ku < size_u() , "Value ku ("+std::to_string(ku)+") should be < size_u ("+std::to_string(size_u())+")");
    ASSERT_CPE(kv >= 0 , "Value kv ("+std::to_string(kv)+") should be >=0 ");
    ASSERT_CPE(kv < size_v() , "Value kv ("+std::to_string(kv)+") should be < size_v ("+std::to_string(size_v())+")");

    int const offset = ku + size_u()*kv;

    ASSERT_CPE(offset < static_cast<int>(speed_data.size()),"Internal error");

    return speed_data[offset];
}

vec3 const& mesh_parametric_cloth::force(int const ku,int const kv) const
{
    ASSERT_CPE(ku >= 0 , "Value ku ("+std::to_string(ku)+") should be >=0 ");
    ASSERT_CPE(ku < size_u() , "Value ku ("+std::to_string(ku)+") should be < size_u ("+std::to_string(size_u())+")");
    ASSERT_CPE(kv >= 0 , "Value kv ("+std::to_string(kv)+") should be >=0 ");
    ASSERT_CPE(kv < size_v() , "Value kv ("+std::to_string(kv)+") should be < size_v ("+std::to_string(size_v())+")");

    int const offset = ku + size_u()*kv;

    ASSERT_CPE(offset < static_cast<int>(force_data.size()),"Internal error");

    return force_data[offset];
}

vec3& mesh_parametric_cloth::force(int const ku,int const kv)
{
    ASSERT_CPE(ku >= 0 , "Value ku ("+std::to_string(ku)+") should be >=0 ");
    ASSERT_CPE(ku < size_u() , "Value ku ("+std::to_string(ku)+") should be < size_u ("+std::to_string(size_u())+")");
    ASSERT_CPE(kv >= 0 , "Value kv ("+std::to_string(kv)+") should be >=0 ");
    ASSERT_CPE(kv < size_v() , "Value kv ("+std::to_string(kv)+") should be < size_v ("+std::to_string(size_v())+")");

    int const offset = ku + size_u()*kv;

    ASSERT_CPE(offset < static_cast<int>(force_data.size()),"Internal error");

    return force_data[offset];
}
vec3 mesh_parametric_cloth::spring_compute_force(int u1, int v1,int u2, int v2, float const L0, float const K)
{
    if (u2<0 || u2>=size_u() || v2<0 || v2>=size_v())
    {
        return vec3(0,0,0);
    }
    else
    {
        /* code */
        // vec3 p0 = vec2(u0,v0);
        vec3 const u0 = vertex(u2,v2)-vertex(u1,v1);
        float const L = norm(u0);
        vec3 const F = K*(L-L0)*u0/L;
        return F;
    }
    
}

void mesh_parametric_cloth::intersect_plane(int ku, int kv){
    if(vertex(ku,kv).z()<-1.100f)
    {
        vertex(ku,kv).z() = -1.100f;
        speed(ku,kv).z() = 0;
        speed(ku,kv) *= 0.9; //on rajoute un coefficient de frottement
    }


}

void mesh_parametric_cloth::intersect_sphere(int ku, int kv, vec3 const& center, float const radius){
    // std::cout<<center[0]<<" "<<center[1]<<" "<<std::endl;
    vec3 cp = vertex(ku,kv) - center;
    if (dot(cp,cp) < radius*radius)
    {
        vec3 n = normalized(cp);
        /* code */
        // std::cout<<"in_cercle"<<std::endl;
        //compute the projection between center and vertex(ku,kv) on the sphere
        vec3 const proj = center + radius * n;
        vec3 vn = dot(speed(ku,kv),n)*n;
        vec3 vt = speed(ku,kv) - vn;
        vertex(ku,kv) = proj;
        speed(ku,kv) = 0.9*vt; //on rajoute un coefficient de frottement
    }
    
}

}
