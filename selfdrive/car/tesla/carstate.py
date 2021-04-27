import copy
from cereal import car,log
from selfdrive.car.tesla.values import DBC, GEAR_MAP, DOORS, BUTTONS, CAR
from selfdrive.car.interfaces import CarStateBase
from opendbc.can.parser import CANParser
from opendbc.can.can_define import CANDefine
from selfdrive.config import Conversions as CV
from selfdrive.car.tesla.CFG_module import read_config_file, load_bool_param

class CarState(CarStateBase):
  def __init__(self, CP):
    super().__init__(CP)
    self.CP = CP
    self.button_states = {button.event_type: False for button in BUTTONS}
    self.can_define = CANDefine(DBC[CP.carFingerprint]['chassis'])

    # Needed by carcontroller
    self.msg_stw_actn_req = None

    self.hands_on_level = 0
    self.das_steeringControl_counter = -1
    self.das_status_counter = -1
    self.steer_warning = None
    self.angle_steers = 0
    self.autopilot_enabled = False
    # 0 = off, 1 = indicate left (stalk down), 2 = indicate right (stalk up)
    self.turn_signal_stalk_state = 0 

    #alca info
    self.prev_alca_pre_engage = False
    self.alca_pre_engage = False
    self.alca_engaged = False
    self.alca_direction = 0
    self.alca_need_engagement = False
    self.alca_done = False

    #to control hazard lighgts
    self.needs_hazard = False

    #path info
    self.laneWidth = 0
    self.lProb = 0.
    self.rProb = 0.
    self.curvC0 = 0.
    self.curvC1 = 0.
    self.curvC2 = 0.
    self.curvC3 = 0.
    self.laneToLeft = 0
    self.laneToRight = 0
    self.lLine = 0
    self.rLine = 0

    self.turn_signal_blinking = False

    #IC integration
    self.DAS_gas_to_resume = 0
    self.DAS_025_steeringOverride = 0 #use for manual steer?
    self.DAS_202_noisyEnvironment = 0 #use for planner error?
    self.DAS_206_apUnavailable = 0 #Ap disabled from CID
    self.DAS_207_lkasUnavailable = 0 #use for manual not in drive?
    self.DAS_208_rackDetected = 0 #use for low battery?
    self.DAS_211_accNoSeatBelt = 0
    self.DAS_216_driverOverriding = 0
    self.DAS_219_lcTempUnavailableSpeed = 0
    self.DAS_220_lcTempUnavailableRoad = 0
    self.DAS_221_lcAborting = 0
    self.DAS_222_accCameraBlind = 0 #we will see what we can use this for
    self.stopSignWarning = 0
    self.stopLightWarning = 0
    self.DAS_canErrors = 0
    self.DAS_notInDrive = 0
    self.DAS_fusedSpeedLimit = 0
    self.baseMapSpeedLimitMPS = 0

    self.cruiseEnabled = False
    self.cruiseDelay = False

    #start config section
    self.forcePedalOverCC = load_bool_param("TinklaEnablePedal",False)
    self.useFollowModeAcc = load_bool_param("TinklaUseFollowACC",False)
    self.autoresumeAcc = load_bool_param("TinklaAutoResumeACC",False)
    self.enableHSO = load_bool_param("TinklaHso",True)
    self.enableHAO = load_bool_param("TinklaHao",False)
    self.enableALC = load_bool_param("TinklaAlc",False)
    self.useTeslaRadar = load_bool_param("TinklaUseTeslaRadar",False)
    self.usesApillarHarness = load_bool_param("TinklaUseAPillarHarness",False)
    self.autoStartAlcaDelay = 2
    self.radarVIN = '"                 "'
    self.radarOffset = 0.0
    self.radarEpasType = 0
    self.radarPosition = 0
    self.hsoNumbPeriod = 1.5
    
    read_config_file(self)
    #end config section

    self.speed_units = "MPH"
    self.tap_direction = 0

    #preAP long
    self.speed_control_enabled = 0
    self.cc_state = 1
    self.adaptive_cruise = 0
    self.apFollowTimeInS = 2.5  # time in seconds to follow
    self.torqueLevel = 0.0
    self.cruise_state = None

    #IC integration
    self.userSpeedLimitOffsetMS = 0

    #pedal interceptor variables
    self.pedal_idx = self.prev_pedal_idx = 0
    self.prev_pedal_interceptor_state = self.pedal_interceptor_state = 0
    self.pedal_interceptor_value = 0.0
    self.pedal_interceptor_value2 = 0.0

  def _convert_to_DAS_fusedSpeedLimit(self, speed_limit_uom, speed_limit_type):
    if speed_limit_uom > 0:
      if speed_limit_type == 0x1E: # Autobahn with no speed limit
        return 0x1F # no speed limit sign
      return int(speed_limit_uom / 5 + 0.5) # sign ID in 5 kph/mph increments (7 shows as 5)
    else:
      if speed_limit_type == 0x1F: # SNA (parking lot, no public road, etc.)
        return 0 # no sign
      return 1 # show 5 kph/mph for unknown limit where we should have one

  def update(self, cp, cp_cam):
    ret = car.CarState.new_message()

    # Vehicle speed
    if (self.CP.carFingerprint == CAR.PREAP_MODELS) and self.usesApillarHarness:
      ret.vEgoRaw = cp_cam.vl["ESP_B"]["ESP_vehicleSpeed"] * CV.KPH_TO_MS
    else:
      ret.vEgoRaw = cp.vl["ESP_B"]["ESP_vehicleSpeed"] * CV.KPH_TO_MS
    ret.vEgo, ret.aEgo = self.update_speed_kf(ret.vEgoRaw)
    ret.standstill = (ret.vEgo < 0.1)

    # Gas pedal
    ret.gas = cp.vl["DI_torque1"]["DI_pedalPos"] / 100.0
    ret.gasPressed = (ret.gas > 0) 
    if self.enableHAO:
      self.DAS_216_driverOverriding = 1 if (ret.gas > 0) else 0
      ret.gas = 0
      ret.gasPressed = False

    # Brake pedal
    ret.brake = 0
    if (self.CP.carFingerprint == CAR.PREAP_MODELS) and self.usesApillarHarness:
      ret.brakePressed = bool(cp_cam.vl["BrakeMessage"]["driverBrakeStatus"] != 1)
    else:
      ret.brakePressed = bool(cp.vl["BrakeMessage"]["driverBrakeStatus"] != 1)
    ret.brakePressed = False #Brake disengages ACC, this just creates issues
    # Steering wheel
    if (self.CP.carFingerprint == CAR.PREAP_MODELS):
      self.hands_on_level = cp_cam.vl["EPAS_sysStatus"]["EPAS_handsOnLevel"]
      self.steer_warning = self.can_define.dv["EPAS_sysStatus"]["EPAS_eacErrorCode"].get(int(cp_cam.vl["EPAS_sysStatus"]["EPAS_eacErrorCode"]), None)
      steer_status = self.can_define.dv["EPAS_sysStatus"]["EPAS_eacStatus"].get(int(cp_cam.vl["EPAS_sysStatus"]["EPAS_eacStatus"]), None)

      ret.steeringAngleDeg = -cp_cam.vl["EPAS_sysStatus"]["EPAS_internalSAS"]
      self.angle_steers = ret.steeringAngleDeg
      ret.steeringTorque = -cp_cam.vl["EPAS_sysStatus"]["EPAS_torsionBarTorque"]
      self.steer_override = abs(cp_cam.vl["EPAS_sysStatus"]["EPAS_handsOnLevel"]) > 0
    else:
      self.hands_on_level = cp.vl["EPAS_sysStatus"]["EPAS_handsOnLevel"]
      self.steer_override = abs(cp.vl["EPAS_sysStatus"]["EPAS_handsOnLevel"]) > 0
      self.steer_warning = self.can_define.dv["EPAS_sysStatus"]["EPAS_eacErrorCode"].get(int(cp.vl["EPAS_sysStatus"]["EPAS_eacErrorCode"]), None)
      steer_status = self.can_define.dv["EPAS_sysStatus"]["EPAS_eacStatus"].get(int(cp.vl["EPAS_sysStatus"]["EPAS_eacStatus"]), None)

      ret.steeringAngleDeg = -cp.vl["EPAS_sysStatus"]["EPAS_internalSAS"]
      self.angle_steers = ret.steeringAngleDeg
      ret.steeringTorque = -cp.vl["EPAS_sysStatus"]["EPAS_torsionBarTorque"]

    ret.steeringRateDeg = -cp.vl["STW_ANGLHP_STAT"]["StW_AnglHP_Spd"] # This is from a different angle sensor, and at different rate

    ret.steeringPressed = (self.hands_on_level > 0)
    ret.steerError = steer_status in ["EAC_FAULT"] #inhibited is not an error, "EAC_INHIBITED"]
    ret.steerWarning = steer_status == "EAC_FAULT" and self.steer_warning in ["EAC_ERROR_MAX_SPEED", "EAC_ERROR_MIN_SPEED", "EAC_ERROR_TMP_FAULT", "SNA"]  # TODO: not sure if this list is complete
    #Trick the alca if autoStartAlcaDelay is set
    if (self.autoStartAlcaDelay > 0) and (self.alca_need_engagement):
      ret.steeringPressed = True
      if self.alca_direction == log.LateralPlan.LaneChangeDirection.left:
        ret.steeringTorque = 1.1
      if self.alca_direction == log.LateralPlan.LaneChangeDirection.right:
        ret.steeringTorque = -1.1

    # Cruise state
    autopilot_status = None
    if (not self.CP.carFingerprint == CAR.PREAP_MODELS):
      autopilot_status = self.can_define.dv["DAS_status"]["DAS_autopilotState"].get(int(cp_cam.vl["DAS_status"]["DAS_autopilotState"]), None)
    autopark_status = None
    eac_status = None
    if (not self.CP.carFingerprint == CAR.PREAP_MODELS):
      eac_status = self.can_define.dv["DAS_pscControl"]["DAS_eacState"].get(int(cp_cam.vl["DAS_pscControl"]["DAS_eacState"]), None)
      autopark_status = self.can_define.dv["DAS_pscControl"]["DAS_pscParkState"].get(int(cp_cam.vl["DAS_pscControl"]["DAS_pscParkState"]), None)
    summon_or_autopark_enabled = (
        eac_status in ["EAC_ACTIVE"] 
        or autopark_status in ["SUMMON", "COMPLETE", "ABORT", "PARALLEL_PULL_OUT_TO_RIGHT", 
        "PARALLEL_PULL_OUT_TO_LEFT", "PARK_RIGHT_CROSS", "PARK_RIGHT_PARALLEL", 
        "PARK_LEFT_CROSS", "PARK_LEFT_PARALLEL"]
    )
    self.cruise_state = cp.vl["DI_state"]["DI_cruiseState"]
    cruise_state = self.can_define.dv["DI_state"]["DI_cruiseState"].get(int(cp.vl["DI_state"]["DI_cruiseState"]), None)
    self.speed_units = self.can_define.dv["DI_state"]["DI_speedUnits"].get(int(cp.vl["DI_state"]["DI_speedUnits"]), None)
    self.v_cruise_actual = cp.vl["DI_state"]["DI_cruiseSet"]
    if self.speed_units == "MPH":
        self.v_cruise_actual = self.v_cruise_actual * CV.MPH_TO_KPH
    self.prev_cruise_buttons = self.cruise_buttons
    self.cruise_buttons = cp.vl["STW_ACTN_RQ"]["SpdCtrlLvr_Stat"]

    acc_enabled = (cruise_state in ["ENABLED", "STANDSTILL", "OVERRIDE", "PRE_FAULT", "PRE_CANCEL"])
    self.autopilot_enabled = (autopilot_status in ["ACTIVE_1", "ACTIVE_2"]) #, "ACTIVE_NAVIGATE_ON_AUTOPILOT"])
    self.cruiseEnabled = acc_enabled and not self.autopilot_enabled and not summon_or_autopark_enabled
    ret.cruiseState.enabled = self.cruiseEnabled and self.cruiseDelay
    if self.speed_units == "KPH":
      ret.cruiseState.speed = cp.vl["DI_state"]["DI_digitalSpeed"] * CV.KPH_TO_MS
    elif self.speed_units == "MPH":
      ret.cruiseState.speed = cp.vl["DI_state"]["DI_digitalSpeed"] * CV.MPH_TO_MS
    ret.cruiseState.available = ((cruise_state == "STANDBY") or ret.cruiseState.enabled)
    ret.cruiseState.standstill = (cruise_state == "STANDSTILL")

    #speed limit
    msu = cp.vl['UI_gpsVehicleSpeed']["UI_mapSpeedLimitUnits"]
    map_speed_uom_to_ms = CV.KPH_TO_MS if msu == 1 else CV.MPH_TO_MS
    map_speed_ms_to_uom = CV.MS_TO_KPH if msu == 1 else CV.MS_TO_MPH
    speed_limit_type = int(cp.vl["UI_driverAssistMapData"]["UI_mapSpeedLimit"])
    rdSignMsg = cp.vl["UI_driverAssistRoadSign"]["UI_roadSign"]
    if rdSignMsg == 3: # ROAD_SIGN_SPEED_LIMIT
      self.baseMapSpeedLimitMPS = cp.vl["UI_driverAssistRoadSign"]["UI_baseMapSpeedLimitMPS"]
      # we round the speed limit in the map's units of measurement to fix noisy data (there are no signs with a limit of 79.2 kph)
      self.baseMapSpeedLimitMPS = int(self.baseMapSpeedLimitMPS * map_speed_ms_to_uom + 0.99) / map_speed_ms_to_uom
    if self.baseMapSpeedLimitMPS > 0 and (speed_limit_type != 0x1F or self.baseMapSpeedLimitMPS <= 5.56):
      self.speed_limit_ms = self.baseMapSpeedLimitMPS # this one is earlier than the actual sign but can also be unreliable, so we ignore it on SNA at higher speeds
    else:
      self.speed_limit_ms = cp.vl['UI_gpsVehicleSpeed']["UI_mppSpeedLimit"] * map_speed_uom_to_ms
    self.DAS_fusedSpeedLimit = self._convert_to_DAS_fusedSpeedLimit(self.speed_limit_ms * map_speed_ms_to_uom, speed_limit_type)

    # Gear
    ret.gearShifter = GEAR_MAP[self.can_define.dv["DI_torque2"]["DI_gear"].get(int(cp.vl["DI_torque2"]["DI_gear"]), "DI_GEAR_INVALID")]
    self.DAS_notInDrive = 0 if ret.gearShifter == car.CarState.GearShifter.drive else 1

    # Buttons
    buttonEvents = []
    for button in BUTTONS:
      state = (cp.vl[button.can_addr][button.can_msg] in button.values)
      if self.button_states[button.event_type] != state:
        event = car.CarState.ButtonEvent.new_message()
        event.type = button.event_type
        event.pressed = state
        buttonEvents.append(event)
      self.button_states[button.event_type] = state
    ret.buttonEvents = buttonEvents

    # Doors
    ret.doorOpen = any([(self.can_define.dv["GTW_carState"][door].get(int(cp.vl["GTW_carState"][door]), "OPEN") == "OPEN") for door in DOORS])

    # Blinkers are used for Comma ALCA
    # we use tap to trigger Comma ALCA so only engage when blinkig but not pressed
    if (cp.vl["STW_ACTN_RQ"]["TurnIndLvr_Stat"] is not None):
      self.turn_signal_stalk_state = (
            0
            if cp.vl["STW_ACTN_RQ"]["TurnIndLvr_Stat"] == 3
            else int(cp.vl["STW_ACTN_RQ"]["TurnIndLvr_Stat"])
        )
    ret.leftBlinker = (cp.vl["GTW_carState"]["BC_indicatorLStatus"] == 1) and (self.turn_signal_stalk_state == 0) and (self.tap_direction == 1)
    ret.rightBlinker = (cp.vl["GTW_carState"]["BC_indicatorRStatus"] == 1) and (self.turn_signal_stalk_state == 0) and (self.tap_direction == 2)
    self.turn_signal_blinking = ret.leftBlinker or ret.rightBlinker

    # Seatbelt
    if (self.CP.carFingerprint == CAR.PREAP_MODELS) and self.usesApillarHarness:
      ret.seatbeltUnlatched = (cp_cam.vl["SDM1"]["SDM_bcklDrivStatus"] != 1)
    elif (self.CP.carFingerprint in [CAR.AP1_MODELX]):
      ret.seatbeltUnlatched = False #TODO: find another way to get this
    else:
      ret.seatbeltUnlatched = (cp.vl["SDM1"]["SDM_bcklDrivStatus"] != 1)

    #Blindspot
    park_right_blindspot = self.can_define.dv["PARK_status2"]["PARK_sdiBlindSpotRight"].get(int(cp.vl["PARK_status2"]["PARK_sdiBlindSpotRight"])) == "WARNING"
    park_left_blindspot = self.can_define.dv["PARK_status2"]["PARK_sdiBlindSpotLeft"].get(int(cp.vl["PARK_status2"]["PARK_sdiBlindSpotLeft"])) == "WARNING"
    if (self.CP.carFingerprint != CAR.PREAP_MODELS):
      das_right_blindspot = self.can_define.dv["DAS_status"]["DAS_blindSpotRearRight"].get(int(cp_cam.vl["DAS_status"]["DAS_blindSpotRearRight"])) in ["WARNING_LEVEL_2","WARNING_LEVEL_1"]
      das_left_blindspot = self.can_define.dv["DAS_status"]["DAS_blindSpotRearLeft"].get(int(cp_cam.vl["DAS_status"]["DAS_blindSpotRearLeft"])) in ["WARNING_LEVEL_2","WARNING_LEVEL_1"]
    else:
      das_right_blindspot = False
      das_left_blindspot = False
    ret.rightBlindspot = (park_right_blindspot or das_right_blindspot)
    ret.leftBlindspot = (park_left_blindspot or das_left_blindspot)
    # Messages needed by carcontroller
    sw_a = copy.copy(cp.vl["STW_ACTN_RQ"])
    if sw_a is not None:
      self.msg_stw_actn_req = sw_a
      
    #Pedal Interceptor
    self.prev_pedal_interceptor_state = self.pedal_interceptor_state
    self.prev_pedal_idx = self.pedal_idx
    if self.CP.carFingerprint in [CAR.PREAP_MODELS]:
      self.pedal_interceptor_state = cp_cam.vl["GAS_SENSOR"]["STATE"]
      self.pedal_interceptor_value = cp_cam.vl["GAS_SENSOR"]["INTERCEPTOR_GAS"]
      self.pedal_interceptor_value2 = cp_cam.vl["GAS_SENSOR"]["INTERCEPTOR_GAS2"]
      self.pedal_idx = cp_cam.vl["GAS_SENSOR"]["IDX"]

    return ret

  @staticmethod
  def get_can_parser(CP):
    usesApillarHarness = load_bool_param("TinklaUseAPillarHarness",False)
    signals = [
      # sig_name, sig_address, default
      ("DI_pedalPos", "DI_torque1", 0),
      ("DI_brakePedal", "DI_torque2", 0),
      ("StW_AnglHP", "STW_ANGLHP_STAT", 0),
      ("StW_AnglHP_Spd", "STW_ANGLHP_STAT", 0),
      ("DI_cruiseState", "DI_state", 0),
      ("DI_digitalSpeed", "DI_state", 0),
      ("DI_speedUnits", "DI_state", 0),
      ("DI_cruiseSet", "DI_state", 0),
      ("DI_gear", "DI_torque2", 0),
      ("DOOR_STATE_FL", "GTW_carState", 1),
      ("DOOR_STATE_FR", "GTW_carState", 1),
      ("DOOR_STATE_RL", "GTW_carState", 1),
      ("DOOR_STATE_RR", "GTW_carState", 1),
      ("DOOR_STATE_FrontTrunk", "GTW_carState", 1),
      ("BOOT_STATE", "GTW_carState", 1),
      ("BC_indicatorLStatus", "GTW_carState", 1),
      ("BC_indicatorRStatus", "GTW_carState", 1),
      # info for speed limit
      ("UI_mapSpeedLimitUnits","UI_gpsVehicleSpeed",0),
      ("UI_mppSpeedLimit","UI_gpsVehicleSpeed",0),
      ("UI_mapSpeedLimit","UI_driverAssistMapData",0),
      ("UI_roadSign","UI_driverAssistRoadSign",0),
      ("UI_baseMapSpeedLimitMPS","UI_driverAssistRoadSign",0),
      # We copy this whole message when spamming cancel
      ("SpdCtrlLvr_Stat", "STW_ACTN_RQ", 0),
      ("VSL_Enbl_Rq", "STW_ACTN_RQ", 0),
      ("SpdCtrlLvrStat_Inv", "STW_ACTN_RQ", 0),
      ("DTR_Dist_Rq", "STW_ACTN_RQ", 0),
      ("TurnIndLvr_Stat", "STW_ACTN_RQ", 0),
      ("HiBmLvr_Stat", "STW_ACTN_RQ", 0),
      ("WprWashSw_Psd", "STW_ACTN_RQ", 0),
      ("WprWash_R_Sw_Posn_V2", "STW_ACTN_RQ", 0),
      ("StW_Lvr_Stat", "STW_ACTN_RQ", 0),
      ("StW_Cond_Flt", "STW_ACTN_RQ", 0),
      ("StW_Cond_Psd", "STW_ACTN_RQ", 0),
      ("HrnSw_Psd", "STW_ACTN_RQ", 0),
      ("StW_Sw00_Psd", "STW_ACTN_RQ", 0),
      ("StW_Sw01_Psd", "STW_ACTN_RQ", 0),
      ("StW_Sw02_Psd", "STW_ACTN_RQ", 0),
      ("StW_Sw03_Psd", "STW_ACTN_RQ", 0),
      ("StW_Sw04_Psd", "STW_ACTN_RQ", 0),
      ("StW_Sw05_Psd", "STW_ACTN_RQ", 0),
      ("StW_Sw06_Psd", "STW_ACTN_RQ", 0),
      ("StW_Sw07_Psd", "STW_ACTN_RQ", 0),
      ("StW_Sw08_Psd", "STW_ACTN_RQ", 0),
      ("StW_Sw09_Psd", "STW_ACTN_RQ", 0),
      ("StW_Sw10_Psd", "STW_ACTN_RQ", 0),
      ("StW_Sw11_Psd", "STW_ACTN_RQ", 0),
      ("StW_Sw12_Psd", "STW_ACTN_RQ", 0),
      ("StW_Sw13_Psd", "STW_ACTN_RQ", 0),
      ("StW_Sw14_Psd", "STW_ACTN_RQ", 0),
      ("StW_Sw15_Psd", "STW_ACTN_RQ", 0),
      ("WprSw6Posn", "STW_ACTN_RQ", 0),
      ("MC_STW_ACTN_RQ", "STW_ACTN_RQ", 0),
      ("CRC_STW_ACTN_RQ", "STW_ACTN_RQ", 0),
      ("PARK_sdiBlindSpotRight","PARK_status2",0),
      ("PARK_sdiBlindSpotLeft","PARK_status2",0),
    ]

    checks = [
      # sig_address, frequency
      ("DI_torque1", 100),
      ("DI_torque2", 100),
      ("STW_ANGLHP_STAT", 100),
      ("DI_state", 10),
      ("STW_ACTN_RQ", 10),
      ("GTW_carState", 10),
    ]

    if not (CP.carFingerprint in [CAR.PREAP_MODELS] and usesApillarHarness):
      signals += [
        ("ESP_vehicleSpeed", "ESP_B", 0),
        ("driverBrakeStatus", "BrakeMessage", 0),
      ]
      checks += [
        ("ESP_B", 50),
        ("BrakeMessage", 50),
      ]
    
    if not (CP.carFingerprint in [CAR.AP1_MODELX]) and not (CP.carFingerprint in [CAR.PREAP_MODELS] and usesApillarHarness):
      signals += [
        ("SDM_bcklDrivStatus", "SDM1", 0),
      ]

      checks += [
        ("SDM1", 10),
      ]

    if CP.carFingerprint in [CAR.AP1_MODELS, CAR.AP2_MODELS, CAR.AP1_MODELX]:
      signals += [
        ("EPAS_handsOnLevel", "EPAS_sysStatus", 0),
        ("EPAS_torsionBarTorque", "EPAS_sysStatus", 0),
        ("EPAS_internalSAS", "EPAS_sysStatus", 0),
        ("EPAS_eacStatus", "EPAS_sysStatus", 1),
        ("EPAS_eacErrorCode", "EPAS_sysStatus", 0),
      ]

      checks += [      
        ("EPAS_sysStatus", 25),
      ]

    return CANParser(DBC[CP.carFingerprint]['chassis'], signals, checks, 0)

  @staticmethod
  def get_cam_can_parser(CP):
    usesApillarHarness = load_bool_param("TinklaUseAPillarHarness",False)
    signals = None
    checks = None
    
    if CP.carFingerprint in [CAR.AP1_MODELS, CAR.AP2_MODELS, CAR.AP1_MODELX]:
      signals = [
        # sig_name, sig_address, default
        #we need the steering control counter
        ("DAS_steeringControlCounter", "DAS_steeringControl", 0),
        # We copy this whole message when we change the status for IC
        ("DAS_autopilotState", "DAS_status", 0),
        ("DAS_blindSpotRearLeft", "DAS_status", 0),
        ("DAS_blindSpotRearRight", "DAS_status", 0),
        ("DAS_fusedSpeedLimit", "DAS_status", 0),
        ("DAS_suppressSpeedWarning", "DAS_status", 0),
        ("DAS_summonObstacle", "DAS_status", 0),
        ("DAS_summonClearedGate", "DAS_status", 0),
        ("DAS_visionOnlySpeedLimit", "DAS_status", 0),
        ("DAS_heaterState", "DAS_status", 0),
        ("DAS_forwardCollisionWarning", "DAS_status", 0),
        ("DAS_autoparkReady", "DAS_status", 0),
        ("DAS_autoParked", "DAS_status", 0),
        ("DAS_autoparkWaitingForBrake", "DAS_status", 0),
        ("DAS_summonFwdLeashReached", "DAS_status", 0),
        ("DAS_summonRvsLeashReached", "DAS_status", 0),
        ("DAS_sideCollisionAvoid", "DAS_status", 0),
        ("DAS_sideCollisionWarning", "DAS_status", 0),
        ("DAS_sideCollisionInhibit", "DAS_status", 0),
        ("DAS_csaState", "DAS_status2", 0),
        ("DAS_laneDepartureWarning", "DAS_status", 0),
        ("DAS_fleetSpeedState", "DAS_status", 0),
        ("DAS_autopilotHandsOnState", "DAS_status", 0),
        ("DAS_autoLaneChangeState", "DAS_status", 0),
        ("DAS_summonAvailable", "DAS_status", 0),
        ("DAS_statusCounter", "DAS_status", 0),
        ("DAS_statusChecksum", "DAS_status", 0),
        ("DAS_headlightRequest", "DAS_bodyControls", 0),
        ("DAS_hazardLightRequest", "DAS_bodyControls", 0),
        ("DAS_wiperSpeed", "DAS_bodyControls", 0),
        ("DAS_turnIndicatorRequest", "DAS_bodyControls", 0),
        ("DAS_highLowBeamDecision", "DAS_bodyControls", 0),
        ("DAS_highLowBeamOffReason", "DAS_bodyControls", 0),
        ("DAS_turnIndicatorRequestReason", "DAS_bodyControls", 0),
        ("DAS_bodyControlsCounter", "DAS_bodyControls", 0),
        ("DAS_bodyControlsChecksum", "DAS_bodyControls", 0),
        ("DAS_accSpeedLimit", "DAS_status2", 0),
        ("DAS_pmmObstacleSeverity", "DAS_status2", 0),
        ("DAS_pmmLoggingRequest", "DAS_status2", 0),
        ("DAS_activationFailureStatus", "DAS_status2", 0),
        ("DAS_pmmUltrasonicsFaultReason", "DAS_status2", 0),
        ("DAS_pmmRadarFaultReason", "DAS_status2", 0),
        ("DAS_pmmSysFaultReason", "DAS_status2", 0),
        ("DAS_pmmCameraFaultReason", "DAS_status2", 0),
        ("DAS_ACC_report", "DAS_status2", 0),
        ("DAS_lssState", "DAS_status", 0),
        ("DAS_radarTelemetry", "DAS_status2", 0),
        ("DAS_robState", "DAS_status2", 0),
        ("DAS_driverInteractionLevel", "DAS_status2", 0),
        ("DAS_ppOffsetDesiredRamp", "DAS_status2", 0),
        ("DAS_longCollisionWarning", "DAS_status2", 0),
        ("DAS_status2Counter", "DAS_status2", 0),
        ("DAS_status2Checksum", "DAS_status2", 0),
        ("DAS_eacState", "DAS_pscControl", 0),
        ("DAS_pscParkState", "DAS_pscControl", 0),
      ]

      checks = [
        # sig_address, frequency
        ("DAS_status", 2),
      ]

    if CP.carFingerprint in [CAR.PREAP_MODELS]:
      signals = [
        ("EPAS_handsOnLevel", "EPAS_sysStatus", 0),
        ("EPAS_torsionBarTorque", "EPAS_sysStatus", 0),
        ("EPAS_internalSAS", "EPAS_sysStatus", 0),
        ("EPAS_eacStatus", "EPAS_sysStatus", 1),
        ("EPAS_eacErrorCode", "EPAS_sysStatus", 0),
        ("INTERCEPTOR_GAS", "GAS_SENSOR", 0),
        ("INTERCEPTOR_GAS2", "GAS_SENSOR", 0),
        ("STATE", "GAS_SENSOR", 0),
        ("IDX", "GAS_SENSOR", 0),
      ]

      checks = [      
        ("EPAS_sysStatus", 25),
      ]

      if usesApillarHarness:
        signals += [
          ("ESP_vehicleSpeed", "ESP_B", 0),
          ("driverBrakeStatus", "BrakeMessage", 0),
          ("SDM_bcklDrivStatus", "SDM1", 0),
        ]
        checks += [
          ("ESP_B", 50),
          ("BrakeMessage", 50),
          ("SDM1", 10),
        ]

    return CANParser(DBC[CP.carFingerprint]['chassis'], signals, checks, 2)
