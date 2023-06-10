#pragma once
#ifndef _FDTD3D_BLEND_ELASTIC_MODEL_MPI_HOST_HPP
#define _FDTD3D_BLEND_ELASTIC_MODEL_MPI_HOST_HPP
#include "fdtd3d_blend_elastic_model_mpi_bones.hpp"
namespace jarvis
{
    inline void fdtd3d_blend_elastic_model_mpi_module<SimulateType::pure_forward>::link()
    {
        glb_model_p = new elasticGeoModel<domain::global>();
        sub_model_p = new elasticGeoModel<domain::local>();
        sub_wave_p = new elasticWave();
        sub_survey_p = new SubGeometry();
        seis_record_p = new SeisRecord();
        sub_pml_p = new mpicuCPML();
        sub_phy_p = new mpicuPhy<SimulateType::pure_forward>();

        sub_pml_p->link(sub_model_p, seis_record_p, sub_wave_p);
        sub_phy_p->link(sub_model_p, seis_record_p, sub_wave_p, sub_survey_p);
        seis_record_p->link(sub_model_p, sub_wave_p, sub_survey_p);
    }

    inline void fdtd3d_blend_elastic_model_mpi_module<SimulateType::pure_forward>::initialize()
    {
        glb_model_p->initialize();
        sub_model_p->initialize(glb_model_p);
        glb_model_p->clear();
        sub_survey_p->initialize(arg_list_p->shot_path, arg_list_p->recv_path, sub_model_p->gridphy);
        seis_record_p->initialize();
        seis_record_p->alloc_host_storage();
        sub_wave_p->initialize_for_org(sub_model_p->gridext);
        sub_pml_p->initialize();
        sub_pml_p->link_geomodel_para(&sub_model_p->lambda, &sub_model_p->mu, &sub_model_p->rho);
        sub_phy_p->initialize();
        sub_phy_p->link_geomodel_para(&sub_model_p->lambda, &sub_model_p->mu, &sub_model_p->rho);
        sub_phy_p->check_is_stable();
        sub_model_p->host_clear();
    }
    inline void fdtd3d_blend_elastic_model_mpi_module<SimulateType::pure_forward>::reset_for_next_shot()
    {
        sub_pml_p->set_zero();
        sub_phy_p->set_zero();
        sub_wave_p->set_zero();
    }
    inline void fdtd3d_blend_elastic_model_mpi_module<SimulateType::pure_forward>::mpi_sub_forward(int _i_shot)
    {
        if (arg_list_p->mpi_rank == 0)
        {
            std::cout << "[forward_modeling_pure]:" << std::endl;
        }
        tic(0);
        for (int it = 0; it < seis_record_p->ntime; it++)
        {
            sub_phy_p->cuda_update_and_mpi_exchange_edge_vel();
            sub_pml_p->cuda_update_and_mpi_exchange_edge_vel();
            sub_phy_p->cuda_update_inside_vel();
            jarvis_mpi_cuda_stream_p->mpi_wait_batch_id(sub_pml_p->mpi_multi_halo_vel.operator_id);
            jarvis_mpi_cuda_stream_p->mpi_wait_batch_id(sub_phy_p->mpi_multi_halo_vel.operator_id);

            sub_phy_p->cuda_add_stress_source(it, _i_shot);
            sub_phy_p->cuda_update_and_mpi_exchange_edge_stress();
            sub_pml_p->cuda_update_and_mpi_exchange_edge_stress();
            sub_phy_p->cuda_update_inside_stress();
            seis_record_p->cuda_get_seis_data(it);
            jarvis_mpi_cuda_stream_p->mpi_wait_batch_id(sub_pml_p->mpi_multi_halo_stress.operator_id);
            jarvis_mpi_cuda_stream_p->mpi_wait_batch_id(sub_phy_p->mpi_multi_halo_stress.operator_id);
            // if (it == 1)
            // {
            //     sub_wave_p->vx.cu_save(arg_list_p->output_path + "fd_vx_1_" + to_string(jarvis_mpi_cuda_stream_p->mpi_frame.mpi_rank), SaveFormat::binary_fld);
            // }
            // if (it == 300)
            // {
            //     sub_wave_p->vx.cu_save(arg_list_p->output_path + "fd_vx_300_" + to_string(jarvis_mpi_cuda_stream_p->mpi_frame.mpi_rank), SaveFormat::binary_fld);
            // }
            // if (it == seis_record_p->ntime - 2)
            // {
            //     std::cout << sizeof(Frame) << std::endl;
            //     sub_wave_p->vx.cu_save(arg_list_p->output_path + "fd_vx_" + to_string(jarvis_mpi_cuda_stream_p->mpi_frame.mpi_rank), SaveFormat::binary_fld);
            // }
        }
        for (int i = 0; i < arg_list_p->mpi_size; i++)
        {
            if (arg_list_p->mpi_rank == i)
            {
                toc(0, "forward modeling pure time:mpi_rank:: " + to_string(arg_list_p->mpi_rank) + ":");
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
    }
    //
    //
    //
    //
    //
    inline void fdtd3d_blend_elastic_model_mpi_module<SimulateType::rtm_forward>::link(elasticGeoModel<domain::local> *_sub_model_p, SubGeometry *_sub_survey_p)
    {
        sub_model_p = _sub_model_p;
        sub_survey_p = _sub_survey_p;
        sub_wave_p = new elasticWave();
        seis_record_p = new SeisRecord();
        sub_pml_p = new mpicuCPML();
        sub_phy_p = new mpicuPhy<SimulateType::rtm_forward>();
        sub_wrc_p = new elasticWaveRecon();
        sub_pml_p->link(sub_model_p, seis_record_p, sub_wave_p);
        sub_phy_p->link(sub_model_p, seis_record_p, sub_wave_p, sub_survey_p);
        sub_wrc_p->link(sub_model_p, seis_record_p, sub_wave_p);
        seis_record_p->link(sub_model_p, sub_wave_p, sub_survey_p);
    }

    inline void fdtd3d_blend_elastic_model_mpi_module<SimulateType::rtm_forward>::initialize()
    {
        seis_record_p->initialize();
        sub_wave_p->initialize_for_rtm(sub_model_p->gridext);
        sub_pml_p->initialize();
        sub_phy_p->initialize();
        sub_wrc_p->initialize();
        sub_pml_p->link_geomodel_para(&sub_model_p->lambda, &sub_model_p->mu, &sub_model_p->rho);
        sub_phy_p->link_geomodel_para(&sub_model_p->lambda, &sub_model_p->mu, &sub_model_p->rho);
    }
    inline void fdtd3d_blend_elastic_model_mpi_module<SimulateType::rtm_forward>::reset_for_next_shot()
    {
        sub_pml_p->set_zero();
        sub_phy_p->set_zero();
        sub_wave_p->set_zero();
    }
    inline void fdtd3d_blend_elastic_model_mpi_module<SimulateType::rtm_forward>::mpi_sub_forward_for_restruct(int _i_shot)
    {
        cuda_tic(0);
        if (arg_list_p->mpi_rank == 0)
            std::cout << "mpi_rank" << arg_list_p->mpi_rank << ":\n\033[42;37m[forward_for_waverecon]\033[0m\ntime loop start..." << std::endl;
        for (int it = 0; it <= seis_record_p->ntime - 1; it++)
        {
            sub_phy_p->cuda_update_and_mpi_exchange_edge_vel();
            sub_pml_p->cuda_update_and_mpi_exchange_edge_vel();
            sub_phy_p->cuda_update_inside_vel();
            jarvis_mpi_cuda_stream_p->mpi_wait_batch_id(sub_pml_p->mpi_multi_halo_vel.operator_id);
            jarvis_mpi_cuda_stream_p->mpi_wait_batch_id(sub_phy_p->mpi_multi_halo_vel.operator_id);

            sub_phy_p->cuda_add_stress_source(it, _i_shot);
            sub_phy_p->cuda_update_and_mpi_exchange_edge_stress();
            sub_pml_p->cuda_update_and_mpi_exchange_edge_stress();
            sub_wrc_p->extract_boundary_wave(it);
            sub_wrc_p->copy_into_host_storage(it);
            sub_phy_p->cuda_update_inside_stress();
            jarvis_mpi_cuda_stream_p->mpi_wait_batch_id(sub_pml_p->mpi_multi_halo_stress.operator_id);
            jarvis_mpi_cuda_stream_p->mpi_wait_batch_id(sub_phy_p->mpi_multi_halo_stress.operator_id);
        }
        for (int i = 0; i < arg_list_p->mpi_size; i++)
        {
            if (arg_list_p->mpi_rank == i)
            {
                cuda_toc(0, "forward for waverecon time::mpi_rank:: " + to_string(arg_list_p->mpi_rank) + ":");
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
    }
    //
    //
    //
    //
    //
    inline void fdtd3d_blend_elastic_model_mpi_module<SimulateType::rtm_backward>::link(fdtd3d_for_restruct_M *_fdtd3d_for_for_restruct_p)
    {
        sub_model_p = _fdtd3d_for_for_restruct_p->sub_model_p;
        sub_wave_p = _fdtd3d_for_for_restruct_p->sub_wave_p;
        sub_survey_p = _fdtd3d_for_for_restruct_p->sub_survey_p;
        sub_wrc_p = _fdtd3d_for_for_restruct_p->sub_wrc_p;
        seis_record_p = _fdtd3d_for_for_restruct_p->seis_record_p;
        sub_phy_p = new mpicuPhy<SimulateType::rtm_backward>();
        sub_phy_p->link(sub_model_p, seis_record_p, sub_wave_p, sub_survey_p);
    }
    inline void fdtd3d_blend_elastic_model_mpi_module<SimulateType::rtm_backward>::initialize()
    {
        sub_phy_p->initialize();
        sub_phy_p->link_geomodel_para(&sub_model_p->lambda, &sub_model_p->mu, &sub_model_p->rho);
    }
    inline void fdtd3d_blend_elastic_model_mpi_module<SimulateType::rtm_backward>::reset_for_next_shot()
    {
        sub_phy_p->set_zero();
        sub_wrc_p->copy_stress_from_host_storage(seis_record_p->ntime - 2);
        sub_wrc_p->copy_vel_from_host_storage(seis_record_p->ntime - 2);
    }
    //
    //
    //
    //
    //
    inline void fdtd3d_blend_elastic_model_mpi_module<SimulateType::rtm_reverse>::link(elasticGeoModel<domain::local> *_sub_model_p, SubGeometry *_sub_survey_p)
    {
        sub_model_p = _sub_model_p;
        sub_survey_p = _sub_survey_p;
        sub_wave_p = new elasticWave();
        seis_record_p = new SeisRecord();
        sub_pml_p = new mpicuCPML();
        sub_phy_p = new mpicuPhy<SimulateType::rtm_reverse>();
        //
        sub_pml_p->link(sub_model_p, seis_record_p, sub_wave_p);
        sub_phy_p->link(sub_model_p, seis_record_p, sub_wave_p, sub_survey_p);
        seis_record_p->link(sub_model_p, sub_wave_p, sub_survey_p);
    }
    inline void fdtd3d_blend_elastic_model_mpi_module<SimulateType::rtm_reverse>::initialize()
    {
        seis_record_p->initialize();
        sub_wave_p->initialize_for_rtm(sub_model_p->gridext);
        sub_pml_p->initialize();
        sub_pml_p->link_geomodel_para(&sub_model_p->lambda, &sub_model_p->mu, &sub_model_p->rho);
        sub_phy_p->initialize();
        sub_phy_p->link_geomodel_para(&sub_model_p->lambda, &sub_model_p->mu, &sub_model_p->rho);
    }
    inline void fdtd3d_blend_elastic_model_mpi_module<SimulateType::rtm_reverse>::reset_for_next_shot()
    {
        sub_pml_p->set_zero();
        sub_phy_p->set_zero();
        sub_wave_p->set_zero();
    }
}
#endif