#pragma once
#ifndef _FDTD3D_GEOMODEL_BONES_HPP
#define _FDTD3D_GEOMODEL_BONES_HPP
#include "../../../../.jarvis_1_util/.jarvis_0_gms/.jarvis_device_gms_header_out.h"
#include "fdtd3d_base_2_wavefield_meat.hpp"
namespace jarvis
{
    enum class domain
    {
        global,
        local
    };
    struct Margin
    {
        int top_margin = 0, bottom_margin = 0;
        int left_margin = 0, right_margin = 0;
        int front_margin = 0, back_margin = 0;
    };
    template <domain domain>
    class elasticGeoModel
    {
    };
    template <>
    class elasticGeoModel<domain::global> : public elastic_fdtd3d_module_Base
    {

    public:
        using host = jarvis::vector<float, MemType::paged>;
        using device = jarvis::vector<float, MemType::device>;
        string model_path;
        Frame gridext;
        Frame gridphy;
        Margin margin;
        //*gms model data
        Field<float, MemType::paged> phy_vp;
        Field<float, MemType::paged> phy_vs;
        Field<float, MemType::paged> phy_rho;

        Field<float, MemType::paged_device> vp;
        Field<float, MemType::paged_device> vs;
        //*org model data
        Field<float, MemType::paged_device> lambda;
        Field<float, MemType::paged_device> mu;
        Field<float, MemType::paged_device> rho;
        //*smooth model data
        Field<float, MemType::paged_device> lambda_smo;
        Field<float, MemType::paged_device> mu_smo;
        Field<float, MemType::paged_device> rho_smo;
        //*inversion model data
        Field<float, MemType::paged_device> lambda_inv;
        Field<float, MemType::paged_device> mu_inv;
        Field<float, MemType::paged_device> rho_inv;
        //
        void initialize();
        void cuda_smooth_model(int x_smo_size, int y_smo_size, int z_smo_size);
        void clear();
    };

    template <>
    class elasticGeoModel<domain::local> : public elastic_fdtd3d_module_Base
    {
    public:
        float glb_vp_min;
        float glb_vs_min;
        float glb_vp_max;
        float glb_vs_max;
        Frame gridext;
        Frame gridphy;
        Margin margin;
        //*org model data
        Field<float, MemType::paged_device> vp;
        Field<float, MemType::paged_device> lambda;
        Field<float, MemType::paged_device> mu;
        Field<float, MemType::paged_device> rho;
        //*smooth model data
        Field<float, MemType::paged_device> lambda_smo;
        Field<float, MemType::paged_device> mu_smo;
        Field<float, MemType::paged_device> rho_smo;
        //*inversion model data
        Field<float, MemType::paged_device> lambda_inv;
        Field<float, MemType::paged_device> mu_inv;
        Field<float, MemType::paged_device> rho_inv;
        //
        void initialize(elasticGeoModel<domain::global> *_glb_model_p);
        void adjust_model_para(float adjust);
        void host_clear();

    private:
        void set_gridext(elasticGeoModel<domain::global> &_glb_model);
        void set_gridphy(elasticGeoModel<domain::global> &_glb_model);
        void set_model_para(const Frame &_glb_gridext,
                            Field<float, MemType::paged_device> &_glb_model_gms_para,
                            Field<float, MemType::paged_device> &_sub_model);
    };
    //
    //
#define set_model_padding(model_pad)              \
    int &top_margin = model_pad.top_margin;       \
    int &bottom_margin = model_pad.bottom_margin; \
    int &left_margin = model_pad.left_margin;     \
    int &right_margin = model_pad.right_margin;   \
    int &front_margin = model_pad.front_margin;   \
    int &back_margin = model_pad.back_margin;
}
#endif