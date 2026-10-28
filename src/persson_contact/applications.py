"""Reserved application interfaces; no leakage or friction solver yet."""


def run_leakage(contact_response, parameters):
    raise NotImplementedError('泄漏求解器尚未实现，当前只预留接口。')


def run_mixed_friction(contact_response, parameters):
    raise NotImplementedError('混合摩擦求解器尚未实现，当前只预留接口。')
