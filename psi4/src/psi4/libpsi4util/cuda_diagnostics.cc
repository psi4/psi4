/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2026 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#ifdef PSI4_HAS_CUDA_RUNTIME

#include "cuda_diagnostics.h"

#include <sstream>

namespace psi {
namespace cuda_diagnostics {

// Some CUDA status codes mean "this GPU or this node is broken/misprovisioned",
// not "Psi4 did something wrong". Reported as a bare "Fatal Error"
// they send users hunting through inputs, keywords, and build flags for a
// defect that does not exist -- an uncorrectable ECC fault, for instance, is
// latched by the driver and fails every subsequent context/stream/allocation on
// that device until it is physically reset, so it can even be inherited from an
// unrelated job that ran earlier on the same card. For those codes, name the
// culprit and give the exact command that confirms it.
//
// Returns an empty string for codes that are not known environmental faults.
std::string cuda_unhealthy_environment_hint(cudaError_t err) {
    switch (err) {
        case cudaErrorECCUncorrectable:
            return
                "  The GPU reported an uncorrectable ECC (double-bit) memory error -- a hardware fault in\n"
                "  the device's own memory. The driver latches this state and then fails every context,\n"
                "  stream, and allocation on that GPU until it is physically reset, so the fault may well\n"
                "  have been triggered by an earlier, unrelated job on the same card.\n"
                "\n"
                "  Confirm with:\n"
                "      nvidia-smi -q -d ECC,ROW_REMAPPER\n"
                "      dmesg | grep -i xid          # look for Xid 48, 63, 92, 94, 95, 171\n"
                "  A non-zero 'DRAM Uncorrectable' count, or 'Remapped Rows ... Pending: Yes', confirms it.\n"
                "\n"
                "  Recovery needs 'nvidia-smi -r' (root, with the GPU idle) or a node reboot. Under a batch\n"
                "  scheduler the faulty device often keeps being handed out: ask an administrator to reset\n"
                "  or drain the node, and resubmit onto different hardware.";
        case cudaErrorDevicesUnavailable:
            return
                "  Every visible GPU is busy or is refusing new contexts. Common causes: the device is in\n"
                "  'Exclusive_Process' compute mode and already held by another process, MIG is enabled but\n"
                "  no instance was assigned, or the scheduler granted a GPU that another job still occupies.\n"
                "\n"
                "  Confirm with:\n"
                "      nvidia-smi --query-gpu=index,compute_mode,mig.mode.current --format=csv\n"
                "      nvidia-smi --query-compute-apps=pid,gpu_uuid,used_memory --format=csv\n"
                "  Also check that CUDA_VISIBLE_DEVICES matches what the scheduler actually granted.";
        case cudaErrorSystemNotReady:
            return
                "  The CUDA system stack has not finished initializing. Typically the NVIDIA fabric manager\n"
                "  is not running, or has not converged, on an NVLink/NVSwitch system.\n"
                "\n"
                "  Confirm with:\n"
                "      systemctl status nvidia-fabricmanager\n"
                "      nvidia-smi -q | grep -i -A2 Fabric\n"
                "  This is a node provisioning problem for an administrator to resolve.";
        case cudaErrorInsufficientDriver:
        case cudaErrorSystemDriverMismatch:
            return
                "  The installed NVIDIA kernel driver is older than, or mismatched with, the CUDA runtime\n"
                "  that Psi4 was built against. This usually means the compute node's driver differs\n"
                "  from the build node's, or a driver upgrade landed without a reboot.\n"
                "\n"
                "  Confirm with:\n"
                "      nvidia-smi --query-gpu=driver_version --format=csv\n"
                "      cat /proc/driver/nvidia/version\n"
                "  Rebuild against the deployed CUDA version, or ask an administrator to align the driver.";
        case cudaErrorNoDevice:
            return
                "  The driver reports no CUDA-capable device at all. Either no GPU was allocated to this\n"
                "  job, CUDA_VISIBLE_DEVICES is empty or names a nonexistent index, or the NVIDIA kernel\n"
                "  module is not loaded on this node.\n"
                "\n"
                "  Confirm with:\n"
                "      echo \"$CUDA_VISIBLE_DEVICES\"; nvidia-smi -L; lsmod | grep nvidia\n"
                "  Under a batch scheduler, verify the job actually requested a GPU (e.g. Slurm --gres=gpu:1).";
        default:
            return std::string();
    }
}

std::string cuda_failure_message(const std::string& what, cudaError_t err, const std::string& device_desc,
                                 const std::string& component) {
    std::ostringstream msg;
    msg << what << ": " << cudaGetErrorString(err);

    const std::string hint = cuda_unhealthy_environment_hint(err);
    if (hint.empty()) return msg.str();

    const std::string blameless = component.empty() ? "Psi4" : "Psi4 or " + component;
    msg << "\n\n"
        << "  This is an unhealthy GPU or node, not a defect in " << blameless << ". No change to the input\n"
        << "  file, keywords, basis set, or Psi4 build will clear it.\n"
        << "\n"
        << "  CUDA status:  " << cudaGetErrorName(err) << " (" << static_cast<int>(err) << ")\n";
    if (!device_desc.empty()) {
        msg << "  Device:       " << device_desc << "\n";
    }
    msg << "\n" << hint << "\n";
    return msg.str();
}

std::string describe_cuda_device(int device_id, const cudaDeviceProp& props) {
    std::ostringstream desc;
    desc << "index " << device_id << " (" << props.name;
    char bus_id[64] = {0};
    if (cudaDeviceGetPCIBusId(bus_id, static_cast<int>(sizeof(bus_id)), device_id) == cudaSuccess && bus_id[0] != '\0') {
        desc << ", PCI " << bus_id;
    }
    desc << ")";
    return desc.str();
}

}  // namespace cuda_diagnostics
}  // namespace psi

#endif  // PSI4_HAS_CUDA_RUNTIME
