#define ESMF_ERR_RETURN(rc) if (ESMF_LogFoundError( \
        rcToCheck=rc, \
        msg=ESMF_LOGERR_PASSTHRU, \
        line=__LINE__, \
        file=__FILE__) \
    ) return

module mpas_nuopc_atm
  !> MPAS NUOPC Cap for Atmosphere

  use ESMF
  use NUOPC
  use NUOPC_Model, &
    modelSS => SetServices

  ! MPAS module required
  use mpas_derived_types, only: core_type, domain_type, block_type, &
       mpas_pool_type, mpas_time_type
  use mpas_kind_types, only: rkind, r8kind, strkind

  implicit none

  private
  ! MPAS variables needed across routines
  ! - mpas_init and mpas_finalize arguments
  type (core_type), pointer :: corelist => null()
  type (domain_type), pointer :: domain => null()
  ! - mpas atm_core_run_{start, advance} variables
  type (block_type), pointer :: block_ptr
  real (kind=rkind), pointer :: dt
  character (len=strkind), pointer :: config_restart_timestamp_name
  real (kind=r8kind) :: diag_start_time, diag_stop_time
  real (kind=r8kind) :: input_start_time, input_stop_time
  real (kind=r8kind) :: output_start_time, output_stop_time
  type (mpas_pool_type), pointer :: state, diag, diag_physics, mesh
  type (mpas_pool_type), pointer :: tend, tend_physics
  logical, pointer :: config_apply_lbcs
  type (MPAS_Time_Type) :: currTime
  character(len=strkind) :: timeStamp
  integer :: itimestep

  character(len=ESMF_MAXSTR), parameter :: file = __FILE__

  public SetVM, SetServices

  !-----------------------------------------------------------------------------
  contains
  !-----------------------------------------------------------------------------

  subroutine SetServices(model, rc)
    !> Register model entry points:
    !>   Advertise: advertise import and export fields
    !>   Realize: realize connected fields
    !>   SetClock: initialize model clock
    !>   DataInitialize: initialize data in import and export states
    !>   Advance: advance model by a single time step
    !>   Finalize: finalize model and cleanup memory

    ! arguments
    type(ESMF_GridComp)  :: model
    integer, intent(out) :: rc

    rc = ESMF_SUCCESS

    ! derive from NUOPC_Model
    call NUOPC_CompDerive(model, modelSS, rc=rc)
    ESMF_ERR_RETURN(rc)

    ! specialize model entry points
    call NUOPC_CompSpecialize(model, specLabel=label_Advertise, &
      specRoutine=Advertise, rc=rc)
    ESMF_ERR_RETURN(rc)
    call NUOPC_CompSpecialize(model, specLabel=label_RealizeProvided, &
      specRoutine=Realize, rc=rc)
    ESMF_ERR_RETURN(rc)
    call NUOPC_CompSpecialize(model, specLabel=label_SetClock, &
      specRoutine=SetClock, rc=rc)
    ESMF_ERR_RETURN(rc)
    call NUOPC_CompSpecialize(model, specLabel=label_DataInitialize, &
      specRoutine=DataInitialize, rc=rc)
    ESMF_ERR_RETURN(rc)
    call NUOPC_CompSpecialize(model, specLabel=label_Advance, &
      specRoutine=Advance, rc=rc)
    ESMF_ERR_RETURN(rc)
    call NUOPC_CompSpecialize(model, specLabel=label_Finalize, &
      specRoutine=Finalize, rc=rc)
    ESMF_ERR_RETURN(rc)

  end subroutine SetServices

  !-----------------------------------------------------------------------------

  subroutine Advertise(model, rc)
    !> Advertise available export fields and desired import fields

    ! arguments
    type(ESMF_GridComp)  :: model
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_VM) :: vm
    type(ESMF_State) :: importState, exportState
    integer :: EXTERNAL_COMM_WORLD

    rc = ESMF_SUCCESS

    call NUOPC_ModelGet(model, importState=importState, &
      exportState=exportState, rc=rc)
    ESMF_ERR_RETURN(rc)

    call ESMF_GridCompGet(model, vm=vm, rc=rc)
    ESMF_ERR_RETURN(rc)

    call ESMF_VMGet(vm, &
      mpiCommunicator=EXTERNAL_COMM_WORLD, rc=rc)
    ESMF_ERR_RETURN(rc)

    call mpas_init(corelist, domain, EXTERNAL_COMM_WORLD)
  end subroutine Advertise

  !-----------------------------------------------------------------------------

  subroutine Realize(model, rc)
    !> Check field connections and realize connected fields

    ! arguments
    type(ESMF_GridComp)  :: model
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_State) :: importState, exportState

    rc = ESMF_SUCCESS

    call NUOPC_ModelGet(model, importState=importState, &
      exportState=exportState, rc=rc)
    ESMF_ERR_RETURN(rc)

    call ESMF_LogWrite(logmsgFlag=ESMF_LOGMSG_ERROR, &
      msg="MPAS NUOPC ATM - Realize has not been implemented", &
      line=__LINE__, &
      file=__FILE__)

  end subroutine Realize

  !-----------------------------------------------------------------------------

  subroutine SetClock(model, rc)
    !> Adjust model clock and time step during initialization
    use atm_core, only: atm_core_run_start

    ! arguments
    type(ESMF_GridComp)  :: model
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Clock) :: clock
    type(ESMF_TimeInterval) :: timeStep
    type(ESMF_Time) :: startTime, currentTime, stopTime
    type(ESMF_State) :: importState, exportState
    logical, pointer :: config_do_restart
    integer :: ierr
    integer(ESMF_KIND_I4) :: dt_i, dt_sec
    character(len=32) :: dateString
    character(len=64) :: msg
    rc = ESMF_SUCCESS

    call NUOPC_ModelGet(model, modelClock=clock, &
      importState=importState, &
      exportState=exportState, rc=rc)
    ESMF_ERR_RETURN(rc)

    ! initialize dt
    call mpas_pool_get_config(domain%blocklist%configs, 'config_dt', dt)
    dt_i = int(dt, kind=ESMF_KIND_I4)
    call ESMF_TimeIntervalSet(timestep, s=dt_i, rc=rc) ! MPAS dt in seconds
    ESMF_ERR_RETURN(rc)
    call ESMF_ClockSet(clock, timeStep=timestep, rc=rc)
    ESMF_ERR_RETURN(rc)

    call NUOPC_CompSetClock(model, clock, rc=rc)
    ESMF_ERR_RETURN(rc)

    call NUOPC_ModelGet(model, modelClock=clock, rc=rc)
    ESMF_ERR_RETURN(rc)
    call ESMF_ClockGet(clock, timestep=timestep, rc=rc)
    call ESMF_TimeIntervalGet(timestep, s=dt_sec, rc=rc)
    write(msg, '(A,I10)') 'ESMF timestep: ', dt_sec
    call ESMF_LogWrite(msg, ESMF_LOGMSG_INFO, rc=rc)

    call ESMF_ClockGet(clock, currTime=currentTime, rc=rc)
    call ESMF_TimeGet(currentTime, timeString=dateString, rc=rc)
    msg = "ESMF current time: " // trim(dateString)
    call ESMF_LogWrite(msg, ESMF_LOGMSG_INFO, rc=rc)

    call ESMF_ClockGet(clock, startTime=startTime, rc=rc)
    call ESMF_TimeGet(startTime, timeString=dateString, rc=rc)
    msg = "ESMF start time: " // trim(dateString)
    call ESMF_LogWrite(msg, ESMF_LOGMSG_INFO, rc=rc)

    call ESMF_ClockGet(clock, stopTime=stopTime, rc=rc)
    call ESMF_TimeGet(stopTime, timeString=dateString, rc=rc)
    msg = "ESMF stop time: " // trim(dateString)
    call ESMF_LogWrite(msg, ESMF_LOGMSG_INFO, rc=rc)

    ierr = atm_core_run_start(domain, block_ptr, dt, config_do_restart, &
         config_restart_timestamp_name, diag_start_time, diag_stop_time, &
         state, diag, diag_physics, mesh, input_start_time, &
         input_stop_time, output_start_time, output_stop_time, &
         config_apply_lbcs, currTime, timestamp, itimestep)

  end subroutine SetClock

  !-----------------------------------------------------------------------------

  subroutine DataInitialize(model, rc)
    !> Initialize data in import and export states

    ! arguments
    type(ESMF_GridComp)  :: model
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Clock) :: clock
    type(ESMF_State) :: importState, exportState

    rc = ESMF_SUCCESS

    call NUOPC_ModelGet(model, modelClock=clock, &
      importState=importState, &
      exportState=exportState, rc=rc)
    ESMF_ERR_RETURN(rc)

    call NUOPC_CompAttributeSet(model, &
      name="InitializeDataComplete", value="true", rc=rc)
    ESMF_ERR_RETURN(rc)

    call ESMF_LogWrite(logmsgFlag=ESMF_LOGMSG_ERROR, &
      msg="MPAS NUOPC ATM - DataInitialize has not been implemented", &
      line=__LINE__, &
      file=__FILE__)

  end subroutine DataInitialize

  !-----------------------------------------------------------------------------

  subroutine Advance(model, rc)
    !> Advance model by a single time step
    use atm_core, only: atm_core_run_advance

    ! arguments
    type(ESMF_GridComp)  :: model
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Clock) :: clock
    type(ESMF_State) :: importState, exportState
    integer :: ierr, stream_dir
    character(len=strkind) :: input_stream, read_time
    real (kind=r8kind) :: integ_start_time, integ_stop_time

    rc = ESMF_SUCCESS

    call NUOPC_ModelGet(model, modelClock=clock, &
      importState=importState, &
      exportState=exportState, rc=rc)
    ESMF_ERR_RETURN(rc)

    ierr = atm_core_run_advance(domain, timestamp, block_ptr, &
         config_apply_lbcs, input_start_time, &
         input_stop_time, output_start_time, output_stop_time, &
         input_stream, read_time, stream_dir, &
         integ_start_time, integ_stop_time, &
         diag_start_time, diag_stop_time, &
         dt, itimestep, state, mesh, diag, diag_physics, &
         tend, tend_physics, config_restart_timestamp_name)

    call ESMF_LogWrite(logmsgFlag=ESMF_LOGMSG_ERROR, &
      msg="MPAS NUOPC ATM - Advance has not been implemented", &
      line=__LINE__, &
      file=__FILE__)

  end subroutine Advance

  !-----------------------------------------------------------------------------

  subroutine Finalize(model, rc)
    !> Finalize model and cleanup memory allocations
    use mpas_subdriver, only: mpas_finalize

    ! arguments
    type(ESMF_GridComp)  :: model
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_State)            :: importState, exportState

    rc = ESMF_SUCCESS

    call NUOPC_ModelGet(model, importState=importState, &
      exportState=exportState, rc=rc)
    ESMF_ERR_RETURN(rc)

    call mpas_finalize(corelist, domain)

  end subroutine Finalize

  !-----------------------------------------------------------------------------

end module mpas_nuopc_atm
