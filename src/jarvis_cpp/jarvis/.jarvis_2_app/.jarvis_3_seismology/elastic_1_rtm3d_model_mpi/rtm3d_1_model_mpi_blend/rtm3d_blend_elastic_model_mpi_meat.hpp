#ifndef _RTM3D_BLEND_ELASTIC_MODEL_MPI_MEAT
#define _RTM3D_BLEND_ELASTIC_MODEL_MPI_MEAT
#include "rtm3d_blend_elastic_model_mpi_bones.hpp"
namespace jarvis
{
    inline void rtm3d_blend_model_mpi_module::link(ArgList *_arg_list_p, jarvis_global_mpicu_Stream_t *_jarvis_mpi_cuda_stream_p)
    {
        arg_list_p = _arg_list_p;
        jarvis_mpi_cuda_stream_p = _jarvis_mpi_cuda_stream_p;
        //
        fdtd_seisdata.link();
        fdtd_rtm_forward.link(fdtd_seisdata.sub_model_p, fdtd_seisdata.sub_survey_p);
        fdtd_rtm_backward.link(&fdtd_rtm_forward);
        fdtd_rtm_reverse.link(fdtd_seisdata.sub_model_p, fdtd_seisdata.sub_survey_p);
        image.link(jarvis_mpi_cuda_stream_p, &fdtd_rtm_backward, &fdtd_rtm_reverse);
    }

    inline void rtm3d_blend_model_mpi_module::mpi_sub_rtm3d()
    {
        fdtd_seisdata.initialize();
        cuda_tic(0);
        for (int ishot = 0; ishot < fdtd_seisdata.sub_survey_p->shot.n_elem; ishot++)
        {
            if (arg_list_p->mpi_rank == 0)
            {
                printf("\033[44;37mshot: %d\033[0m\n", ishot + 1);
            }
            fdtd_seisdata.reset_for_next_shot();
            fdtd_seisdata.mpi_sub_forward(ishot);
            fdtd_seisdata.seis_record_p->copy_seis_into_host(ishot);
            if (fdtd_seisdata.sub_survey_p->is_have_recv)
            {
                fdtd_seisdata.seis_record_p->vx_seis.cu_save(arg_list_p->output_path + "vx_seis_" + to_string(arg_list_p->mpi_rank) + "_" + to_string(ishot), SaveFormat::ascii_txt);
                fdtd_seisdata.seis_record_p->vy_seis.cu_save(arg_list_p->output_path + "vy_seis_" + to_string(arg_list_p->mpi_rank) + "_" + to_string(ishot), SaveFormat::ascii_txt);
                fdtd_seisdata.seis_record_p->vz_seis.cu_save(arg_list_p->output_path + "vz_seis_" + to_string(arg_list_p->mpi_rank) + "_" + to_string(ishot), SaveFormat::ascii_txt);
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
        cuda_toc(0, " ");
        fdtd_seisdata.sub_wave_p->clear();
        fdtd_seisdata.sub_pml_p->clear();
        fdtd_seisdata.sub_phy_p->clear();
        jarvis_mpi_cuda_stream_p->clear();

        fdtd_rtm_forward.initialize();
        fdtd_rtm_backward.initialize();
        fdtd_rtm_reverse.initialize();
        image.initialize();
        // tic(0);
        for (int ishot = 0; ishot < fdtd_seisdata.sub_survey_p->shot.n_elem; ishot++)
        {
            if (arg_list_p->mpi_rank == 0)
            {
                printf("\033[44;37mshot: %d\033[0m\n", ishot + 1);
            }
            fdtd_seisdata.seis_record_p->copy_seis_from_host(ishot);
            fdtd_rtm_forward.reset_for_next_shot();
            fdtd_rtm_forward.mpi_sub_forward_for_restruct(ishot);
            //
            fdtd_rtm_backward.reset_for_next_shot();
            fdtd_rtm_reverse.reset_for_next_shot();
            cuda_tic(0);
            //* TIME REVERSE
            for (int it = fdtd_seisdata.seis_record_p->ntime - 1; it >= 0; it--)
            {
                //*UPDATE 1
                fdtd_rtm_backward.sub_phy_p->cuda_subtract_stress_source(it, ishot);
                fdtd_rtm_reverse.sub_phy_p->cuda_rtm_reverse_add_vel_source(it, fdtd_seisdata.seis_record_p);
                //
                fdtd_rtm_backward.sub_phy_p->cuda_update_and_mpi_exchange_edge_stress();
                //
                fdtd_rtm_reverse.sub_phy_p->cuda_update_and_mpi_exchange_edge_vel();
                fdtd_rtm_reverse.sub_pml_p->cuda_update_and_mpi_exchange_edge_vel();
                //
                fdtd_rtm_backward.sub_phy_p->cuda_update_inside_stress();
                fdtd_rtm_backward.sub_wrc_p->inject_boundary_stress(it - 1);
                fdtd_rtm_backward.sub_wrc_p->copy_stress_from_host_storage(it - 2);
                //
                fdtd_rtm_reverse.sub_phy_p->cuda_update_inside_vel();
                //
                jarvis_mpi_cuda_stream_p->mpi_wait_batch_id(fdtd_rtm_reverse.sub_pml_p->mpi_multi_halo_vel.operator_id);
                jarvis_mpi_cuda_stream_p->mpi_wait_batch_id(fdtd_rtm_reverse.sub_phy_p->mpi_multi_halo_vel.operator_id);
                jarvis_mpi_cuda_stream_p->mpi_wait_batch_id(fdtd_rtm_backward.sub_phy_p->mpi_multi_halo_stress.operator_id);
                //*UPDATE 2
                fdtd_rtm_backward.sub_phy_p->cuda_update_and_mpi_exchange_edge_vel();
                //
                fdtd_rtm_reverse.sub_phy_p->cuda_update_and_mpi_exchange_edge_stress();
                fdtd_rtm_reverse.sub_pml_p->cuda_update_and_mpi_exchange_edge_stress();
                //
                fdtd_rtm_backward.sub_phy_p->cuda_update_inside_vel();
                fdtd_rtm_backward.sub_wrc_p->inject_boundary_vel(it - 1);
                fdtd_rtm_backward.sub_wrc_p->copy_vel_from_host_storage(it - 2);
                //
                fdtd_rtm_reverse.sub_phy_p->cuda_update_inside_stress();
                //
                image.cuda_cal_correlation_from_wave();
                jarvis_mpi_cuda_stream_p->mpi_wait_batch_id(fdtd_rtm_reverse.sub_pml_p->mpi_multi_halo_stress.operator_id);
                jarvis_mpi_cuda_stream_p->mpi_wait_batch_id(fdtd_rtm_reverse.sub_phy_p->mpi_multi_halo_stress.operator_id);
                jarvis_mpi_cuda_stream_p->mpi_wait_batch_id(fdtd_rtm_backward.sub_phy_p->mpi_multi_halo_vel.operator_id);

                // if (it == 1)
                // {
                //     fdtd_rtm_backward.sub_wave_p->vx.cu_save(arg_list_p->output_path + "wrc_vx_1_" + to_string(jarvis_mpi_cuda_stream_p->mpi_frame.mpi_rank), SaveFormat::binary_fld);
                // }
                // if (it == 301)
                // {
                //     fdtd_rtm_backward.sub_wave_p->vx.cu_save(arg_list_p->output_path + "wrc_vx_300_" + to_string(jarvis_mpi_cuda_stream_p->mpi_frame.mpi_rank), SaveFormat::binary_fld);
                // }
                // if (it == 500)
                // {
                //     fdtd_rtm_backward.sub_wave_p->vx.cu_save(arg_list_p->output_path + "wrc_vx_500_" + to_string(jarvis_mpi_cuda_stream_p->mpi_frame.mpi_rank), SaveFormat::binary_fld);
                // }
            }
            for (int i = 0; i < arg_list_p->mpi_size; i++)
            {
                if (arg_list_p->mpi_rank == i)
                {
                    cuda_toc(0, "rtm of one shot time::mpi_rank:: " + to_string(arg_list_p->mpi_rank) + ":");
                }
                MPI_Barrier(MPI_COMM_WORLD);
            }
        }
        std::cout << arg_list_p->mpi_rank << std::endl;
        for (int i = 0; i < arg_list_p->mpi_size; i++)
        {
            if (arg_list_p->mpi_rank == i)
            {
                image.image_pp.cu_save(arg_list_p->output_path + "image_pp_" + to_string(arg_list_p->mpi_rank), SaveFormat::binary_raw);
                image.Ipp_laplace.cu_save(arg_list_p->output_path + "image_pp_laplace_" + to_string(arg_list_p->mpi_rank), SaveFormat::binary_raw);
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
    }
}
#endif
