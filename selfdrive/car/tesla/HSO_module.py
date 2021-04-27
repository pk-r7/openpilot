# human steer override module

class HSOController:
    def __init__(self, carcontroller):
        self.human_control = False
        self.frame_humanSteered = 0

    def update_stat(self, CC, CS, enabled, actuators, frame):
        human_control = False

        if CS.enableHSO and enabled:
            # if steering but not by ALCA
            if CS.steer_override > 0:
                self.frame_humanSteered = frame
            elif (frame - self.frame_humanSteered < 50) and (CS.turn_signal_stalk_state > 0):  
                # stalk locked, update frame
                self.frame_humanSteered = frame
            elif (frame - self.frame_humanSteered < 50):  
                # Need more human testing of handoff timing
                # Find steering difference between visiond model and human (no need to do every frame if we run out of CPU):
                apply_steer = int(actuators.steeringAngleDeg)
                angle_diff = abs(apply_steer - CS.angle_steers)
                if angle_diff > 15.0:
                    self.frame_humanSteered = frame
            if frame - self.frame_humanSteered < 50:
                human_control = True

        self.human_control = human_control
        human_control = human_control and enabled
        hands_on_fault = (CS.hands_on_level >= 2 and not human_control)
        lkas_enabled = enabled and (not hands_on_fault)
        return human_control, hands_on_fault, lkas_enabled
