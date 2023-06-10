#ifndef _RTM3D_BLEND_ELASTIC_MODEL_MPI_BONES
#define _RTM3D_BLEND_ELASTIC_MODEL_MPI_BONES
#include ".rtm3d_blend_elastic_model_mpi_header_in.h"
namespace jarvis
{
	class rtm3d_blend_model_mpi_module : public elastic_fdtd3d_module_Base
	{
	public:
		fdtd3d_blend_elastic_model_mpi_module<SimulateType::pure_forward> fdtd_seisdata;
		fdtd3d_blend_elastic_model_mpi_module<SimulateType::rtm_forward> fdtd_rtm_forward;
		fdtd3d_blend_elastic_model_mpi_module<SimulateType::rtm_backward> fdtd_rtm_backward;
		fdtd3d_blend_elastic_model_mpi_module<SimulateType::rtm_reverse> fdtd_rtm_reverse;
		//
		MigImage image;
		// ADCIGs adcig;
		void link(ArgList *_arg_list_p, jarvis_global_mpicu_Stream_t *_jarvis_mpi_cuda_stream_p);
		//
		void mpi_sub_rtm3d();
	};

	void device_mpi_sub_rtm3d(rtm3d_blend_model_mpi_module *rtm_mpi);
}
#endif