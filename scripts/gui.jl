"""
    SLAM.jl GUI - Scattering Lagrange Asymptotic Matching

A graphical interface for nuclear scattering calculations using
Lagrange mesh methods with Baye's exact differential matrices.

Usage:
    julia gui.jl

Requirements:
    - Blink.jl
    - JSON3
"""

using Blink
using JSON3
using LinearAlgebra
using Printf

# Set BLAS threads for performance
BLAS.set_num_threads(Sys.CPU_THREADS)

# Include SLAM module
include(joinpath(@__DIR__, "..", "src", "SLAM.jl"))
using .SLAM

# Global variables
global calculation_running = false
global current_slam_results = nothing
global current_numerov_results = nothing
global current_wavefunctions = nothing

# HTML/CSS/JavaScript for the UI
const HTML_CONTENT = """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>SLAM Nuclear Scattering Calculator</title>
    <link rel="icon" type="image/png" href="data:image/png;base64,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">
    <style>
        * {
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }

        body {
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, sans-serif;
            background: linear-gradient(135deg, #1a5276 0%, #2e86de 100%);
            min-height: 100vh;
            padding: 20px;
        }

        .container {
            max-width: 1800px;
            margin: 0 auto;
        }

        .header {
            background: linear-gradient(135deg, #0a1628 0%, #1a3a5c 100%);
            padding: 15px 25px;
            border-radius: 15px;
            margin-bottom: 20px;
            box-shadow: 0 8px 32px rgba(0, 0, 0, 0.3);
            display: flex;
            align-items: center;
            justify-content: center;
            gap: 20px;
        }

        .header-logo {
            width: 100px;
            height: 100px;
            flex-shrink: 0;
        }

        .header-text {
            text-align: left;
        }

        .header h1 {
            color: #00d4ff;
            font-size: 2.0em;
            font-weight: 700;
            margin-bottom: 4px;
            text-shadow: 0 0 10px rgba(0, 212, 255, 0.5);
        }

        .header p {
            color: #5dade2;
            font-size: 0.95em;
            letter-spacing: 1px;
        }

        .main-grid {
            display: grid;
            grid-template-columns: 380px 1fr 380px;
            gap: 20px;
            margin-bottom: 20px;
        }

        .card {
            background: rgba(255, 255, 255, 0.95);
            border-radius: 15px;
            padding: 20px;
            box-shadow: 0 8px 32px rgba(0, 0, 0, 0.1);
        }

        .card-title {
            color: #1a5276;
            font-size: 1.3em;
            font-weight: 600;
            margin-bottom: 15px;
            padding-bottom: 10px;
            border-bottom: 2px solid #2e86de;
        }

        .param-section {
            margin-bottom: 20px;
        }

        .section-label {
            color: #2e86de;
            font-weight: 600;
            font-size: 1.0em;
            margin-bottom: 12px;
        }

        .param-row {
            display: flex;
            align-items: center;
            margin-bottom: 10px;
        }

        .param-label {
            flex: 0 0 90px;
            color: #555;
            font-weight: 500;
            font-size: 0.9em;
        }

        .param-value {
            flex: 0 0 60px;
            text-align: center;
            color: #1a5276;
            font-weight: 600;
            font-size: 1.0em;
        }

        input[type="range"] {
            flex: 1;
            margin: 0 8px;
            -webkit-appearance: none;
            height: 6px;
            border-radius: 3px;
            background: linear-gradient(90deg, #1a5276, #2e86de);
            outline: none;
        }

        input[type="range"]::-webkit-slider-thumb {
            -webkit-appearance: none;
            width: 16px;
            height: 16px;
            border-radius: 50%;
            background: #fff;
            cursor: pointer;
            box-shadow: 0 2px 8px rgba(0,0,0,0.2);
        }

        input[type="number"] {
            flex: 1;
            margin: 0 8px;
            padding: 6px 10px;
            border: 2px solid #ddd;
            border-radius: 6px;
            font-size: 0.9em;
            color: #555;
            max-width: 100px;
        }

        input[type="number"]:focus {
            border-color: #2e86de;
            outline: none;
        }

        select {
            width: 100%;
            padding: 8px;
            border: 2px solid #ddd;
            border-radius: 6px;
            font-size: 0.9em;
            color: #555;
            background: white;
            cursor: pointer;
            transition: all 0.3s;
        }

        select:hover {
            border-color: #2e86de;
        }

        .toggle-container {
            display: flex;
            align-items: center;
            justify-content: space-between;
            margin-bottom: 10px;
        }

        .toggle {
            position: relative;
            width: 50px;
            height: 26px;
            background: #ddd;
            border-radius: 13px;
            cursor: pointer;
            transition: all 0.3s;
        }

        .toggle.active {
            background: linear-gradient(90deg, #1a5276, #2e86de);
        }

        .toggle-slider {
            position: absolute;
            top: 3px;
            left: 3px;
            width: 20px;
            height: 20px;
            background: white;
            border-radius: 50%;
            transition: all 0.3s;
        }

        .toggle.active .toggle-slider {
            left: 27px;
        }

        .run-button {
            width: 100%;
            padding: 12px;
            background: linear-gradient(135deg, #1a5276, #2e86de);
            color: white;
            border: none;
            border-radius: 8px;
            font-size: 1.1em;
            font-weight: 600;
            cursor: pointer;
            transition: all 0.3s;
            box-shadow: 0 4px 15px rgba(26, 82, 118, 0.4);
            margin-bottom: 8px;
        }

        .run-button:hover {
            transform: translateY(-2px);
            box-shadow: 0 6px 20px rgba(26, 82, 118, 0.6);
        }

        .run-button:active {
            transform: translateY(0);
        }

        .run-button:disabled {
            background: #ccc;
            cursor: not-allowed;
            transform: none;
        }

        .status {
            text-align: center;
            padding: 8px;
            margin-top: 8px;
            border-radius: 6px;
            font-weight: 500;
            font-size: 0.9em;
        }

        .status.ready { background: #e8f5e9; color: #2e7d32; }
        .status.running { background: #fff3e0; color: #e65100; }
        .status.success { background: #e8f5e9; color: #2e7d32; }
        .status.error { background: #ffebee; color: #c62828; }

        .output-console {
            height: 450px;
            overflow-y: auto;
            background: #1e1e1e;
            color: #d4d4d4;
            padding: 12px;
            border-radius: 8px;
            font-family: 'Consolas', 'Monaco', monospace;
            font-size: 0.85em;
            line-height: 1.4;
        }

        .output-console::-webkit-scrollbar {
            width: 8px;
        }

        .output-console::-webkit-scrollbar-track {
            background: #2d2d2d;
        }

        .output-console::-webkit-scrollbar-thumb {
            background: #2e86de;
            border-radius: 4px;
        }

        .result-item {
            padding: 12px;
            background: #f5f5f5;
            border-left: 4px solid #2e86de;
            margin-bottom: 8px;
            border-radius: 5px;
        }

        .result-label {
            color: #666;
            font-size: 0.85em;
            margin-bottom: 4px;
        }

        .result-value {
            color: #333;
            font-size: 1.1em;
            font-weight: 600;
            font-family: 'Consolas', 'Monaco', monospace;
        }

        .result-comparison {
            display: grid;
            grid-template-columns: 1fr 1fr;
            gap: 8px;
            padding: 12px;
            background: #f0f8ff;
            border-left: 4px solid #1a5276;
            margin-bottom: 8px;
            border-radius: 5px;
        }

        .result-comparison .method-label {
            font-size: 0.75em;
            color: #666;
            margin-bottom: 2px;
        }

        .result-comparison .method-value {
            font-size: 0.95em;
            font-weight: 600;
            font-family: 'Consolas', 'Monaco', monospace;
        }

        .slam-val { color: #1a5276; }
        .numerov-val { color: #e74c3c; }

        .viz-container {
            background: rgba(255, 255, 255, 0.95);
            border-radius: 15px;
            padding: 20px;
            box-shadow: 0 8px 32px rgba(0, 0, 0, 0.1);
        }

        .plot-area {
            display: grid;
            grid-template-columns: 1fr 1fr;
            gap: 20px;
            margin-top: 15px;
        }

        .plot-canvas {
            width: 100%;
            height: 380px;
            background: #f9f9f9;
            border-radius: 10px;
        }

        .plot-controls {
            display: flex;
            align-items: center;
            gap: 15px;
            margin-bottom: 10px;
        }

        .plot-controls label {
            font-weight: 500;
            color: #555;
        }

        .plot-controls select {
            width: auto;
            min-width: 150px;
        }

        @keyframes spin {
            0% { transform: rotate(0deg); }
            100% { transform: rotate(360deg); }
        }

        .spinner {
            border: 4px solid #f3f3f3;
            border-top: 4px solid #2e86de;
            border-radius: 50%;
            width: 40px;
            height: 40px;
            animation: spin 1s linear infinite;
            margin: 0 auto;
        }

        .tab-container {
            display: flex;
            margin-bottom: 15px;
            border-bottom: 2px solid #eee;
        }

        .tab {
            padding: 8px 16px;
            cursor: pointer;
            border-bottom: 2px solid transparent;
            margin-bottom: -2px;
            color: #666;
            font-weight: 500;
            transition: all 0.3s;
        }

        .tab:hover {
            color: #2e86de;
        }

        .tab.active {
            color: #1a5276;
            border-bottom-color: #2e86de;
        }

        .tab-content {
            display: none;
        }

        .tab-content.active {
            display: block;
        }

        .potential-grid {
            display: grid;
            grid-template-columns: 1fr 1fr 1fr;
            gap: 8px;
        }

        .potential-input {
            display: flex;
            flex-direction: column;
        }

        .potential-input label {
            font-size: 0.75em;
            color: #666;
            margin-bottom: 2px;
        }

        .potential-input input {
            padding: 4px 6px;
            border: 1px solid #ddd;
            border-radius: 4px;
            font-size: 0.85em;
            width: 100%;
        }
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <div class="header-logo">
                <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 200 200" width="100" height="100">
                  <defs>
                    <radialGradient id="nucleusGrad" cx="50%" cy="50%" r="50%">
                      <stop offset="0%" style="stop-color:#5dade2"/>
                      <stop offset="50%" style="stop-color:#2e86de"/>
                      <stop offset="100%" style="stop-color:#1a5276"/>
                    </radialGradient>
                    <radialGradient id="glowGrad" cx="50%" cy="50%" r="50%">
                      <stop offset="0%" style="stop-color:#3498db;stop-opacity:0.8"/>
                      <stop offset="100%" style="stop-color:#3498db;stop-opacity:0"/>
                    </radialGradient>
                    <filter id="glow" x="-100%" y="-100%" width="300%" height="300%">
                      <feGaussianBlur stdDeviation="3" result="coloredBlur"/>
                      <feMerge><feMergeNode in="coloredBlur"/><feMergeNode in="SourceGraphic"/></feMerge>
                    </filter>
                    <filter id="strongGlow" x="-100%" y="-100%" width="300%" height="300%">
                      <feGaussianBlur stdDeviation="4" result="coloredBlur"/>
                      <feMerge><feMergeNode in="coloredBlur"/><feMergeNode in="coloredBlur"/><feMergeNode in="SourceGraphic"/></feMerge>
                    </filter>
                  </defs>
                  <!-- Animated glow -->
                  <circle cx="100" cy="100" r="40" fill="url(#glowGrad)">
                    <animate attributeName="r" values="35;50;35" dur="2s" repeatCount="indefinite"/>
                    <animate attributeName="opacity" values="0.4;0.8;0.4" dur="2s" repeatCount="indefinite"/>
                  </circle>
                  <!-- Rotating partial wave circles -->
                  <g stroke="#1e90ff" fill="none" opacity="0.5">
                    <circle cx="100" cy="100" r="55" stroke-dasharray="8,8" stroke-width="1">
                      <animateTransform attributeName="transform" type="rotate" from="0 100 100" to="360 100 100" dur="15s" repeatCount="indefinite"/>
                    </circle>
                    <circle cx="100" cy="100" r="75" stroke-dasharray="12,6" stroke-width="1">
                      <animateTransform attributeName="transform" type="rotate" from="360 100 100" to="0 100 100" dur="20s" repeatCount="indefinite"/>
                    </circle>
                  </g>
                  <!-- Incoming pulses (hidden by default, shown during calculation) -->
                  <g id="incoming-particles" style="display:none;">
                    <!-- Bright projectile balls -->
                    <circle r="8" fill="#00ff88" filter="url(#strongGlow)">
                      <animate attributeName="cx" values="-10;100" dur="0.8s" repeatCount="indefinite" begin="0s"/>
                      <animate attributeName="cy" values="100;100" dur="0.8s" repeatCount="indefinite" begin="0s"/>
                      <animate attributeName="opacity" values="1;1;0" dur="0.8s" repeatCount="indefinite" begin="0s"/>
                    </circle>
                    <circle r="8" fill="#00ff88" filter="url(#strongGlow)">
                      <animate attributeName="cx" values="-10;100" dur="0.8s" repeatCount="indefinite" begin="0.27s"/>
                      <animate attributeName="cy" values="100;100" dur="0.8s" repeatCount="indefinite" begin="0.27s"/>
                      <animate attributeName="opacity" values="1;1;0" dur="0.8s" repeatCount="indefinite" begin="0.27s"/>
                    </circle>
                    <circle r="8" fill="#00ff88" filter="url(#strongGlow)">
                      <animate attributeName="cx" values="-10;100" dur="0.8s" repeatCount="indefinite" begin="0.54s"/>
                      <animate attributeName="cy" values="100;100" dur="0.8s" repeatCount="indefinite" begin="0.54s"/>
                      <animate attributeName="opacity" values="1;1;0" dur="0.8s" repeatCount="indefinite" begin="0.54s"/>
                    </circle>
                    <!-- Trail lines -->
                    <line stroke="#00ff88" stroke-width="4" filter="url(#glow)">
                      <animate attributeName="x1" values="-30;80" dur="0.8s" repeatCount="indefinite" begin="0s"/>
                      <animate attributeName="x2" values="-10;100" dur="0.8s" repeatCount="indefinite" begin="0s"/>
                      <animate attributeName="y1" values="100;100" dur="0.8s" repeatCount="indefinite" begin="0s"/>
                      <animate attributeName="y2" values="100;100" dur="0.8s" repeatCount="indefinite" begin="0s"/>
                      <animate attributeName="opacity" values="0.8;0.8;0" dur="0.8s" repeatCount="indefinite" begin="0s"/>
                    </line>
                    <line stroke="#00ff88" stroke-width="4" filter="url(#glow)">
                      <animate attributeName="x1" values="-30;80" dur="0.8s" repeatCount="indefinite" begin="0.27s"/>
                      <animate attributeName="x2" values="-10;100" dur="0.8s" repeatCount="indefinite" begin="0.27s"/>
                      <animate attributeName="y1" values="100;100" dur="0.8s" repeatCount="indefinite" begin="0.27s"/>
                      <animate attributeName="y2" values="100;100" dur="0.8s" repeatCount="indefinite" begin="0.27s"/>
                      <animate attributeName="opacity" values="0.8;0.8;0" dur="0.8s" repeatCount="indefinite" begin="0.27s"/>
                    </line>
                    <line stroke="#00ff88" stroke-width="4" filter="url(#glow)">
                      <animate attributeName="x1" values="-30;80" dur="0.8s" repeatCount="indefinite" begin="0.54s"/>
                      <animate attributeName="x2" values="-10;100" dur="0.8s" repeatCount="indefinite" begin="0.54s"/>
                      <animate attributeName="y1" values="100;100" dur="0.8s" repeatCount="indefinite" begin="0.54s"/>
                      <animate attributeName="y2" values="100;100" dur="0.8s" repeatCount="indefinite" begin="0.54s"/>
                      <animate attributeName="opacity" values="0.8;0.8;0" dur="0.8s" repeatCount="indefinite" begin="0.54s"/>
                    </line>
                  </g>
                  <!-- Expanding scatter rings -->
                  <circle cx="100" cy="100" r="30" stroke="#00d4ff" stroke-width="2" fill="none" opacity="0">
                    <animate attributeName="r" values="30;95" dur="1.8s" repeatCount="indefinite"/>
                    <animate attributeName="opacity" values="0.8;0" dur="1.8s" repeatCount="indefinite"/>
                    <animate attributeName="stroke-width" values="2;0.3" dur="1.8s" repeatCount="indefinite"/>
                  </circle>
                  <circle cx="100" cy="100" r="30" stroke="#00d4ff" stroke-width="2" fill="none" opacity="0">
                    <animate attributeName="r" values="30;95" dur="1.8s" repeatCount="indefinite" begin="0.6s"/>
                    <animate attributeName="opacity" values="0.8;0" dur="1.8s" repeatCount="indefinite" begin="0.6s"/>
                    <animate attributeName="stroke-width" values="2;0.3" dur="1.8s" repeatCount="indefinite" begin="0.6s"/>
                  </circle>
                  <circle cx="100" cy="100" r="30" stroke="#00d4ff" stroke-width="2" fill="none" opacity="0">
                    <animate attributeName="r" values="30;95" dur="1.8s" repeatCount="indefinite" begin="1.2s"/>
                    <animate attributeName="opacity" values="0.8;0" dur="1.8s" repeatCount="indefinite" begin="1.2s"/>
                    <animate attributeName="stroke-width" values="2;0.3" dur="1.8s" repeatCount="indefinite" begin="1.2s"/>
                  </circle>
                  <!-- Nucleus -->
                  <circle cx="100" cy="100" r="28" fill="url(#nucleusGrad)" filter="url(#strongGlow)"/>
                  <!-- Nucleons -->
                  <circle cx="92" cy="92" r="7" fill="#e74c3c" opacity="0.9">
                    <animate attributeName="cx" values="92;95;92" dur="2.5s" repeatCount="indefinite"/>
                    <animate attributeName="cy" values="92;94;92" dur="2.8s" repeatCount="indefinite"/>
                  </circle>
                  <circle cx="108" cy="94" r="7" fill="#5dade2" opacity="0.9">
                    <animate attributeName="cx" values="108;105;108" dur="2.7s" repeatCount="indefinite"/>
                    <animate attributeName="cy" values="94;92;94" dur="2.4s" repeatCount="indefinite"/>
                  </circle>
                  <circle cx="90" cy="108" r="7" fill="#5dade2" opacity="0.9">
                    <animate attributeName="cx" values="90;93;90" dur="2.6s" repeatCount="indefinite"/>
                    <animate attributeName="cy" values="108;105;108" dur="2.9s" repeatCount="indefinite"/>
                  </circle>
                  <circle cx="106" cy="106" r="7" fill="#e74c3c" opacity="0.9">
                    <animate attributeName="cx" values="106;103;106" dur="2.8s" repeatCount="indefinite"/>
                    <animate attributeName="cy" values="106;109;106" dur="2.5s" repeatCount="indefinite"/>
                  </circle>
                  <!-- S symbol -->
                  <text x="100" y="107" font-family="Georgia, serif" font-size="24" font-weight="bold" font-style="italic" fill="white" text-anchor="middle" filter="url(#glow)">S</text>
                </svg>
            </div>
            <div class="header-text">
                <h1>SLAM Nuclear Scattering Calculator</h1>
                <p>Scattering Lagrange Asymptotic Matching</p>
            </div>
        </div>

        <div class="main-grid">
            <!-- Left Panel: Parameters -->
            <div class="card">
                <div class="card-title">Parameters</div>

                <div class="tab-container">
                    <div class="tab active" data-tab="scattering">Scattering</div>
                    <div class="tab" data-tab="potential">Potential</div>
                    <div class="tab" data-tab="numerical">Numerical</div>
                </div>

                <!-- Scattering Tab -->
                <div id="scattering-tab" class="tab-content active">
                    <div class="param-section">
                        <div class="section-label">Energy & Angular Momentum</div>
                        <div class="param-row">
                            <span class="param-label">E_lab (MeV):</span>
                            <input type="number" id="e_lab" value="30.0" step="1" min="0.1">
                        </div>
                        <div class="param-row">
                            <span class="param-label">l_max:</span>
                            <input type="range" id="l_max" min="5" max="30" step="1" value="15">
                            <span class="param-value" id="l_max-val">15</span>
                        </div>
                        <div class="param-row">
                            <span class="param-label">Spin:</span>
                            <select id="spin_val">
                                <option value="0">0 (spinless)</option>
                                <option value="0.5" selected>1/2</option>
                                <option value="1">1</option>
                                <option value="1.5">3/2</option>
                            </select>
                        </div>
                    </div>

                    <div class="param-section">
                        <div class="section-label">Target / Projectile</div>
                        <div class="param-row">
                            <span class="param-label">Z_proj:</span>
                            <input type="number" id="z_proj" value="1" step="1" min="0">
                        </div>
                        <div class="param-row">
                            <span class="param-label">A_proj:</span>
                            <input type="number" id="a_proj" value="1" step="1" min="1">
                        </div>
                        <div class="param-row">
                            <span class="param-label">Z_targ:</span>
                            <input type="number" id="z_targ" value="20" step="1" min="0">
                        </div>
                        <div class="param-row">
                            <span class="param-label">A_targ:</span>
                            <input type="number" id="a_targ" value="40" step="1" min="1">
                        </div>
                    </div>
                </div>

                <!-- Potential Tab -->
                <div id="potential-tab" class="tab-content">
                    <div class="param-section">
                        <div class="section-label">Volume (Real)</div>
                        <div class="potential-grid">
                            <div class="potential-input">
                                <label>V_v (MeV)</label>
                                <input type="number" id="V_v" value="53.0" step="0.1">
                            </div>
                            <div class="potential-input">
                                <label>r_v (fm)</label>
                                <input type="number" id="r_v" value="1.25" step="0.01">
                            </div>
                            <div class="potential-input">
                                <label>a_v (fm)</label>
                                <input type="number" id="a_v" value="0.65" step="0.01">
                            </div>
                        </div>
                    </div>

                    <div class="param-section">
                        <div class="section-label">Volume (Imaginary)</div>
                        <div class="potential-grid">
                            <div class="potential-input">
                                <label>W_v (MeV)</label>
                                <input type="number" id="W_v" value="0.0" step="0.1">
                            </div>
                            <div class="potential-input">
                                <label>r_wv (fm)</label>
                                <input type="number" id="r_wv" value="1.25" step="0.01">
                            </div>
                            <div class="potential-input">
                                <label>a_wv (fm)</label>
                                <input type="number" id="a_wv" value="0.65" step="0.01">
                            </div>
                        </div>
                    </div>

                    <div class="param-section">
                        <div class="section-label">Surface (Imaginary)</div>
                        <div class="potential-grid">
                            <div class="potential-input">
                                <label>W_s (MeV)</label>
                                <input type="number" id="W_s" value="10.0" step="0.1">
                            </div>
                            <div class="potential-input">
                                <label>r_ws (fm)</label>
                                <input type="number" id="r_ws" value="1.25" step="0.01">
                            </div>
                            <div class="potential-input">
                                <label>a_ws (fm)</label>
                                <input type="number" id="a_ws" value="0.47" step="0.01">
                            </div>
                        </div>
                    </div>

                    <div class="param-section">
                        <div class="section-label">Spin-Orbit (Real)</div>
                        <div class="potential-grid">
                            <div class="potential-input">
                                <label>V_so (MeV)</label>
                                <input type="number" id="V_so" value="6.0" step="0.1">
                            </div>
                            <div class="potential-input">
                                <label>r_so (fm)</label>
                                <input type="number" id="r_so" value="1.10" step="0.01">
                            </div>
                            <div class="potential-input">
                                <label>a_so (fm)</label>
                                <input type="number" id="a_so" value="0.65" step="0.01">
                            </div>
                        </div>
                    </div>

                    <div class="param-section">
                        <div class="section-label">Coulomb</div>
                        <div class="param-row">
                            <span class="param-label">r_c (fm):</span>
                            <input type="number" id="r_c" value="1.25" step="0.01">
                        </div>
                    </div>
                </div>

                <!-- Numerical Tab -->
                <div id="numerical-tab" class="tab-content">
                    <div class="param-section">
                        <div class="section-label">SLAM Method</div>
                        <div class="param-row">
                            <span class="param-label">N (points):</span>
                            <input type="range" id="N_mesh" min="20" max="100" step="10" value="50">
                            <span class="param-value" id="N_mesh-val">50</span>
                        </div>
                        <div class="param-row">
                            <span class="param-label">R (fm):</span>
                            <input type="number" id="R_max" value="15.0" step="1" min="5">
                        </div>
                    </div>

                    <div class="param-section">
                        <div class="section-label">Numerov Method</div>
                        <div class="param-row">
                            <span class="param-label">h (fm):</span>
                            <input type="number" id="h_step" value="0.05" step="0.01" min="0.01" max="0.2">
                        </div>
                    </div>

                    <div class="param-section">
                        <div class="toggle-container">
                            <span class="section-label">Compare with Numerov</span>
                            <div class="toggle active" id="compare-toggle">
                                <div class="toggle-slider"></div>
                            </div>
                        </div>
                    </div>
                </div>

                <button class="run-button" id="run-btn">Run All Partial Waves</button>
                <div class="status ready" id="status">Ready</div>
            </div>

            <!-- Middle Panel: Output -->
            <div class="card">
                <div class="card-title">Calculation Output</div>
                <div class="output-console" id="output"></div>
            </div>

            <!-- Right Panel: Results -->
            <div class="card">
                <div class="card-title">Results Summary</div>
                <div id="results">
                    <div class="result-item">
                        <div class="result-label">Total Partial Waves</div>
                        <div class="result-value" id="n-waves">--</div>
                    </div>
                    <div class="result-item">
                        <div class="result-label">Wave Number k (fm^-1)</div>
                        <div class="result-value" id="k-val">--</div>
                    </div>
                    <div class="result-item">
                        <div class="result-label">Coulomb Parameter eta</div>
                        <div class="result-value" id="eta-val">--</div>
                    </div>
                    <div class="result-item">
                        <div class="result-label">Total Reaction Cross Section</div>
                        <div class="result-value" id="cross-section">--</div>
                    </div>
                    <div class="result-comparison">
                        <div>
                            <div class="method-label">Max |S| Diff (SLAM vs Numerov)</div>
                            <div class="method-value slam-val" id="max-s-diff">--</div>
                        </div>
                        <div>
                            <div class="method-label">Avg |S| Diff</div>
                            <div class="method-value numerov-val" id="avg-s-diff">--</div>
                        </div>
                    </div>
                </div>
            </div>
        </div>

        <!-- Bottom: Visualization -->
        <div class="viz-container">
            <div class="card-title">Visualization</div>
            <div class="plot-controls">
                <label>Select Partial Wave for Wavefunction:</label>
                <select id="wave-select">
                    <option value="0">l=0, j=0.5</option>
                </select>
            </div>
            <div class="plot-area">
                <div id="plot1" class="plot-canvas" style="display:none; width:100%; height:450px;"></div>
                <div id="plot2" class="plot-canvas" style="display:none; width:100%; height:450px;"></div>
            </div>
        </div>
    </div>

    <script>
        // Wait for Plotly to load
        function waitForPlotly(callback, maxAttempts = 50) {
            var attempts = 0;
            var checkInterval = setInterval(function() {
                attempts++;
                if (typeof Plotly !== 'undefined') {
                    console.log('Plotly loaded successfully after', attempts, 'attempts');
                    clearInterval(checkInterval);
                    callback();
                } else if (attempts >= maxAttempts) {
                    console.error('Failed to load Plotly after', maxAttempts, 'attempts');
                    clearInterval(checkInterval);
                }
            }, 100);
        }

        // Global data storage
        window.wavefunction_data = {};
        window.smatrix_data = {};

        // Store wavefunction data (called from Julia)
        window.storeWaveData = function(key, r_slam, psi_slam_re, psi_slam_im, r_num, psi_num_re, psi_num_im) {
            window.wavefunction_data[key] = {
                r_slam: r_slam,
                psi_slam_re: psi_slam_re,
                psi_slam_im: psi_slam_im,
                r_num: r_num,
                psi_num_re: psi_num_re,
                psi_num_im: psi_num_im
            };
            console.log('Stored wavefunction data for key:', key);
        };

        // Update wave select dropdown
        window.updateWaveSelect = function(options) {
            var select = document.getElementById('wave-select');
            select.innerHTML = '';
            options.forEach(function(opt) {
                var option = document.createElement('option');
                option.value = opt.value;
                option.textContent = opt.label;
                select.appendChild(option);
            });
        };

        // Plot wavefunction for selected partial wave
        window.plotWavefunction = function(l, j, r_slam, psi_slam_re, psi_slam_im, r_num, psi_num_re, psi_num_im) {
            if (typeof Plotly === 'undefined') {
                waitForPlotly(function() {
                    window.plotWavefunction(l, j, r_slam, psi_slam_re, psi_slam_im, r_num, psi_num_re, psi_num_im);
                });
                return;
            }

            // Store data for later use
            var key = l + '_' + j;
            window.wavefunction_data[key] = {
                r_slam: r_slam, psi_slam_re: psi_slam_re, psi_slam_im: psi_slam_im,
                r_num: r_num, psi_num_re: psi_num_re, psi_num_im: psi_num_im
            };

            try {
                var traces = [
                    {
                        x: r_slam, y: psi_slam_re,
                        name: 'Re(ψ) SLAM',
                        type: 'scatter', mode: 'lines',
                        line: {width: 2, color: '#1a5276'}
                    },
                    {
                        x: r_slam, y: psi_slam_im,
                        name: 'Im(ψ) SLAM',
                        type: 'scatter', mode: 'lines',
                        line: {width: 2, color: '#2e86de', dash: 'dash'}
                    },
                    {
                        x: r_num, y: psi_num_re,
                        name: 'Re(ψ) Numerov',
                        type: 'scatter', mode: 'lines',
                        line: {width: 2, color: '#e74c3c'}
                    },
                    {
                        x: r_num, y: psi_num_im,
                        name: 'Im(ψ) Numerov',
                        type: 'scatter', mode: 'lines',
                        line: {width: 2, color: '#c0392b', dash: 'dash'}
                    }
                ];

                var layout = {
                    title: 'Wave Function (l=' + l + ', j=' + j + ')',
                    xaxis: {title: 'r (fm)'},
                    yaxis: {title: 'ψ(r)'},
                    height: 430,
                    legend: {x: 0.65, y: 0.95}
                };

                Plotly.newPlot('plot1', traces, layout);
                document.getElementById('plot1').style.display = 'block';
            } catch(e) {
                console.error('Error in plotWavefunction:', e);
            }
        };

        // Plot |S| comparison for all partial waves
        window.plotSMatrixComparison = function(l_values, j_values, s_slam, s_numerov) {
            if (typeof Plotly === 'undefined') {
                waitForPlotly(function() {
                    window.plotSMatrixComparison(l_values, j_values, s_slam, s_numerov);
                });
                return;
            }

            // Store data
            window.smatrix_data = {l: l_values, j: j_values, s_slam: s_slam, s_numerov: s_numerov};

            try {
                // Create labels like "0,0.5", "1,0.5", "1,1.5", etc.
                var labels = [];
                for (var i = 0; i < l_values.length; i++) {
                    labels.push('(' + l_values[i] + ',' + j_values[i] + ')');
                }

                var traces = [
                    {
                        x: labels, y: s_slam,
                        name: 'SLAM',
                        type: 'scatter', mode: 'lines+markers',
                        line: {width: 2, color: '#1a5276'},
                        marker: {size: 8}
                    },
                    {
                        x: labels, y: s_numerov,
                        name: 'Numerov',
                        type: 'scatter', mode: 'lines+markers',
                        line: {width: 2, color: '#e74c3c', dash: 'dash'},
                        marker: {size: 6}
                    }
                ];

                var layout = {
                    title: '|S| Comparison: SLAM vs Numerov',
                    xaxis: {title: 'Partial Wave (l, j)', tickangle: -45},
                    yaxis: {title: '|S|', range: [0, 1.1]},
                    height: 430,
                    legend: {x: 0.75, y: 0.95}
                };

                Plotly.newPlot('plot2', traces, layout);
                document.getElementById('plot2').style.display = 'block';
            } catch(e) {
                console.error('Error in plotSMatrixComparison:', e);
            }
        };

        // Update wavefunction plot when selection changes
        document.getElementById('wave-select').addEventListener('change', function() {
            var key = this.value;
            console.log('Selected wave:', key);
            console.log('Available data keys:', Object.keys(window.wavefunction_data));
            var data = window.wavefunction_data[key];
            if (data) {
                var parts = key.split('_');
                var l = parseInt(parts[0]);
                var j = parseFloat(parts[1]);
                console.log('Plotting l=', l, 'j=', j);

                // Directly create the plot instead of calling plotWavefunction to ensure update
                var traces = [
                    {
                        x: data.r_slam, y: data.psi_slam_re,
                        name: 'Re(ψ) SLAM',
                        type: 'scatter', mode: 'lines',
                        line: {width: 2, color: '#1a5276'}
                    },
                    {
                        x: data.r_slam, y: data.psi_slam_im,
                        name: 'Im(ψ) SLAM',
                        type: 'scatter', mode: 'lines',
                        line: {width: 2, color: '#2e86de', dash: 'dash'}
                    },
                    {
                        x: data.r_num, y: data.psi_num_re,
                        name: 'Re(ψ) Numerov',
                        type: 'scatter', mode: 'lines',
                        line: {width: 2, color: '#e74c3c'}
                    },
                    {
                        x: data.r_num, y: data.psi_num_im,
                        name: 'Im(ψ) Numerov',
                        type: 'scatter', mode: 'lines',
                        line: {width: 2, color: '#c0392b', dash: 'dash'}
                    }
                ];

                var layout = {
                    title: 'Wave Function (l=' + l + ', j=' + j + ')',
                    xaxis: {title: 'r (fm)'},
                    yaxis: {title: 'ψ(r)'},
                    height: 430,
                    legend: {x: 0.65, y: 0.95}
                };

                Plotly.react('plot1', traces, layout);
            } else {
                console.log('No data found for key:', key);
            }
        });

        // Tab switching
        document.querySelectorAll('.tab').forEach(tab => {
            tab.addEventListener('click', function() {
                document.querySelectorAll('.tab').forEach(t => t.classList.remove('active'));
                document.querySelectorAll('.tab-content').forEach(c => c.classList.remove('active'));

                this.classList.add('active');
                document.getElementById(this.dataset.tab + '-tab').classList.add('active');
            });
        });

        // Update slider values
        const sliders = ['l_max', 'N_mesh'];
        sliders.forEach(id => {
            const slider = document.getElementById(id);
            const display = document.getElementById(id + '-val');
            slider.addEventListener('input', () => {
                display.textContent = slider.value;
            });
        });

        // Toggle switch
        let compareEnabled = true;
        document.getElementById('compare-toggle').addEventListener('click', function() {
            compareEnabled = !compareEnabled;
            this.classList.toggle('active');
        });

        // Output functions
        function appendOutput(text) {
            const output = document.getElementById('output');
            output.innerHTML += text + '<br>';
            output.scrollTop = output.scrollHeight;
        }

        function clearOutput() {
            document.getElementById('output').innerHTML = '';
        }

        function setStatus(text, className) {
            const status = document.getElementById('status');
            status.textContent = text;
            status.className = 'status ' + className;
        }

        // Logo animation control
        function showIncomingParticles() {
            var particles = document.getElementById('incoming-particles');
            if (particles) particles.style.display = 'block';
        }

        function hideIncomingParticles() {
            var particles = document.getElementById('incoming-particles');
            if (particles) particles.style.display = 'none';
        }

        // Export for Julia to call
        window.showIncomingParticles = showIncomingParticles;
        window.hideIncomingParticles = hideIncomingParticles;

        // Run calculation
        document.getElementById('run-btn').addEventListener('click', async function() {
            const btn = this;
            btn.disabled = true;
            clearOutput();
            setStatus('Running calculation...', 'running');

            // Show incoming particles animation
            showIncomingParticles();

            // Clear previous data
            window.wavefunction_data = {};
            window.smatrix_data = {};

            const params = {
                e_lab: parseFloat(document.getElementById('e_lab').value),
                l_max: parseInt(document.getElementById('l_max').value),
                spin: parseFloat(document.getElementById('spin_val').value),
                z_proj: parseFloat(document.getElementById('z_proj').value),
                a_proj: parseFloat(document.getElementById('a_proj').value),
                z_targ: parseFloat(document.getElementById('z_targ').value),
                a_targ: parseFloat(document.getElementById('a_targ').value),
                V_v: parseFloat(document.getElementById('V_v').value),
                r_v: parseFloat(document.getElementById('r_v').value),
                a_v: parseFloat(document.getElementById('a_v').value),
                W_v: parseFloat(document.getElementById('W_v').value),
                r_wv: parseFloat(document.getElementById('r_wv').value),
                a_wv: parseFloat(document.getElementById('a_wv').value),
                W_s: parseFloat(document.getElementById('W_s').value),
                r_ws: parseFloat(document.getElementById('r_ws').value),
                a_ws: parseFloat(document.getElementById('a_ws').value),
                V_so: parseFloat(document.getElementById('V_so').value),
                r_so: parseFloat(document.getElementById('r_so').value),
                a_so: parseFloat(document.getElementById('a_so').value),
                r_c: parseFloat(document.getElementById('r_c').value),
                N: parseInt(document.getElementById('N_mesh').value),
                R: parseFloat(document.getElementById('R_max').value),
                h: parseFloat(document.getElementById('h_step').value),
                compare: compareEnabled
            };

            try {
                await Blink.msg("run_calculation", params);
            } catch(e) {
                appendOutput("<br><br>Error: " + e);
                setStatus('Communication error', 'error');
                btn.disabled = false;
            }
        });
    </script>
</body>
</html>
"""

function run_calculation(params::Dict)
    global win, current_slam_results, current_numerov_results, current_wavefunctions

    try
        global calculation_running = true

        # Extract parameters (ensure correct types)
        E_lab = Float64(params["e_lab"])
        l_max = Int(params["l_max"])
        spin = Float64(params["spin"])
        Z_proj = Float64(params["z_proj"])
        A_proj = Float64(params["a_proj"])
        Z_targ = Float64(params["z_targ"])
        A_targ = Float64(params["a_targ"])
        N = Int(params["N"])
        R = Float64(params["R"])
        h = Float64(params["h"])
        compare = params["compare"]

        # Calculate CM energy
        E_cm = E_lab * A_targ / (A_proj + A_targ)

        # Output function
        function web_output(msg)
            fullmsg = msg * "<br>"
            @js_ win document.getElementById("output").innerHTML += $fullmsg
            @js_ win document.getElementById("output").scrollTop = document.getElementById("output").scrollHeight
        end

        web_output("="^60)
        web_output("    SLAM SCATTERING CALCULATION")
        web_output("="^60)
        web_output("")
        web_output(@sprintf("Projectile: Z=%d, A=%d", Int(Z_proj), Int(A_proj)))
        web_output(@sprintf("Target:     Z=%d, A=%d", Int(Z_targ), Int(A_targ)))
        web_output(@sprintf("E_lab = %.2f MeV  ->  E_cm = %.4f MeV", E_lab, E_cm))
        web_output(@sprintf("Spin = %.1f, l_max = %d", spin, l_max))
        web_output("")

        # Create optical potential
        pot = OpticalPotential(
            V_v=Float64(params["V_v"]), r_v=Float64(params["r_v"]), a_v=Float64(params["a_v"]),
            W_v=Float64(params["W_v"]), r_wv=Float64(params["r_wv"]), a_wv=Float64(params["a_wv"]),
            W_s=Float64(params["W_s"]), r_ws=Float64(params["r_ws"]), a_ws=Float64(params["a_ws"]),
            V_so=Float64(params["V_so"]), r_so=Float64(params["r_so"]), a_so=Float64(params["a_so"]),
            r_c=Float64(params["r_c"]),
            Z_proj=Z_proj, A_proj=A_proj,
            Z_targ=Z_targ, A_targ=A_targ
        )

        web_output("Optical Potential Parameters:")
        web_output(@sprintf("  Volume:  V_v=%.2f MeV, r_v=%.2f fm, a_v=%.2f fm", params["V_v"], params["r_v"], params["a_v"]))
        web_output(@sprintf("  Surface: W_s=%.2f MeV, r_ws=%.2f fm, a_ws=%.2f fm", params["W_s"], params["r_ws"], params["a_ws"]))
        web_output(@sprintf("  Spin-orbit: V_so=%.2f MeV, r_so=%.2f fm, a_so=%.2f fm", params["V_so"], params["r_so"], params["a_so"]))
        web_output(@sprintf("  Coulomb: r_c=%.2f fm", params["r_c"]))
        web_output("")
        web_output(@sprintf("Numerical: N=%d, R=%.1f fm, h=%.3f fm", N, R, h))
        web_output("")
        web_output("-"^60)

        # Generate j values based on spin
        function get_j_values(l::Int, s::Float64)
            if s == 0.0
                return [Float64(l)]
            elseif l == 0
                return [s]
            else
                return [abs(l - s), l + s]
            end
        end

        # Storage for results
        slam_results = []
        numerov_results = []
        wavefunctions = Dict()

        l_values = Float64[]
        j_values = Float64[]
        s_slam = Float64[]
        s_numerov = Float64[]
        wave_options = []

        web_output("Computing all partial waves...")
        web_output("")
        web_output("<span style='font-family: monospace;'>┌─────┬──────┬──────────┬──────────┬──────────┐</span>")
        web_output("<span style='font-family: monospace;'>│  l  │  j   │ |S| SLAM │|S| Numrv │   Diff   │</span>")
        web_output("<span style='font-family: monospace;'>├─────┼──────┼──────────┼──────────┼──────────┤</span>")

        first_result = nothing

        for l in 0:l_max
            for j in get_j_values(l, spin)
                # SLAM calculation
                prob = ScatteringProblem(pot, E_cm, l; j=j, N=N, R=R)
                slam_res = solve_scattering(prob)
                push!(slam_results, slam_res)

                if first_result === nothing
                    first_result = slam_res
                end

                # Numerov calculation if enabled
                if compare
                    num_res = solve_numerov(prob; h=h)
                    push!(numerov_results, num_res)

                    s_diff = abs(abs(slam_res.S_matrix) - abs(num_res.S_matrix))
                    web_output(@sprintf("<span style='font-family: monospace;'>│ %2d  │ %4.1f │ %.6f │ %.6f │ %.2e │</span>",
                                       l, j, abs(slam_res.S_matrix), abs(num_res.S_matrix), s_diff))

                    push!(s_numerov, abs(num_res.S_matrix))

                    # Store wavefunction for plotting
                    wf = solve_wavefunction(prob)
                    wavefunctions["$(l)_$(j)"] = (wf, num_res)
                else
                    web_output(@sprintf("<span style='font-family: monospace;'>│ %2d  │ %4.1f │ %.6f │    --    │    --    │</span>", l, j, abs(slam_res.S_matrix)))
                    push!(s_numerov, 0.0)
                end

                push!(l_values, Float64(l))
                push!(j_values, j)
                push!(s_slam, abs(slam_res.S_matrix))
                push!(wave_options, Dict("value" => "$(l)_$(j)", "label" => "l=$l, j=$j"))
            end
        end

        current_slam_results = slam_results
        current_numerov_results = numerov_results
        current_wavefunctions = wavefunctions

        # Close table
        web_output("<span style='font-family: monospace;'>└─────┴──────┴──────────┴──────────┴──────────┘</span>")
        web_output("")

        # Total reaction cross section
        σ_R = 0.0
        k = first_result.k
        η = first_result.η
        for res in slam_results
            σ_R += (2*res.j + 1) * (1 - abs2(res.S_matrix))
        end
        σ_R *= π / k^2 * 10.0  # Convert to mb

        web_output(@sprintf("Total Reaction Cross Section: %.2f mb", σ_R))

        # Update results display
        n_waves_text = string(length(slam_results))
        k_text = @sprintf("%.6f", k)
        eta_text = @sprintf("%.6f", η)
        σ_R_text = @sprintf("%.2f mb", σ_R)

        @js_ win document.getElementById("n-waves").textContent = $n_waves_text
        @js_ win document.getElementById("k-val").textContent = $k_text
        @js_ win document.getElementById("eta-val").textContent = $eta_text
        @js_ win document.getElementById("cross-section").textContent = $σ_R_text

        if compare && length(numerov_results) > 0
            s_diffs = abs.(s_slam .- s_numerov)
            max_diff = maximum(s_diffs)
            avg_diff = mean(s_diffs)

            max_diff_text = @sprintf("%.2e", max_diff)
            avg_diff_text = @sprintf("%.2e", avg_diff)

            @js_ win document.getElementById("max-s-diff").textContent = $max_diff_text
            @js_ win document.getElementById("avg-s-diff").textContent = $avg_diff_text
        end

        # Update wave select dropdown
        @js_ win window.updateWaveSelect($wave_options)

        # Plot |S| comparison
        web_output("")
        web_output("Generating plots...")

        @js_ win window.plotSMatrixComparison($l_values, $j_values, $s_slam, $s_numerov)

        # Send all wavefunction data to JavaScript for selection
        if compare && length(wavefunctions) > 0
            web_output("Transferring wavefunction data...")

            for (key, (wf, num_res)) in wavefunctions
                # SLAM wavefunction
                r_slam = collect(wf.r)
                psi_slam_re = real.(wf.ψ_total)
                psi_slam_im = imag.(wf.ψ_total)

                # Numerov wavefunction
                r_num = collect(num_res.r)
                psi_num_re = real.(num_res.ψ_normalized)
                psi_num_im = imag.(num_res.ψ_normalized)

                # Store data in JavaScript using function call
                @js_ win window.storeWaveData($key, $r_slam, $psi_slam_re, $psi_slam_im,
                                               $r_num, $psi_num_re, $psi_num_im)
            end

            # Plot first wavefunction (l=0)
            first_key = "0_$(spin == 0.0 ? 0.0 : spin)"
            if haskey(wavefunctions, first_key)
                wf, num_res = wavefunctions[first_key]

                r_slam = collect(wf.r)
                psi_slam_re = real.(wf.ψ_total)
                psi_slam_im = imag.(wf.ψ_total)

                r_num = collect(num_res.r)
                psi_num_re = real.(num_res.ψ_normalized)
                psi_num_im = imag.(num_res.ψ_normalized)

                l_first = 0
                j_first = spin == 0.0 ? 0.0 : spin

                @js_ win window.plotWavefunction($l_first, $j_first, $r_slam, $psi_slam_re, $psi_slam_im,
                                                  $r_num, $psi_num_re, $psi_num_im)
            end
        end

        web_output("")
        web_output("="^60)
        web_output("    CALCULATION COMPLETE")
        web_output("="^60)

        # Hide incoming particles animation
        @js_ win window.hideIncomingParticles()

        @js_ win document.getElementById("status").textContent = "Calculation complete!"
        @js_ win document.getElementById("status").className = "status success"
        @js_ win document.getElementById("run-btn").disabled = false

    catch e
        # Hide incoming particles animation on error
        @js_ win window.hideIncomingParticles()

        @js_ win document.getElementById("status").textContent = "Error occurred"
        @js_ win document.getElementById("status").className = "status error"

        error_msg = sprint(showerror, e)
        full_error = "<br>ERROR:<br>" * replace(error_msg, "\n" => "<br>")
        @js_ win document.getElementById("output").innerHTML += $full_error

        println("Error: ", e)
        println(stacktrace(catch_backtrace()))

        @js_ win document.getElementById("run-btn").disabled = false
    finally
        global calculation_running = false
    end
end

# Helper function for mean
function mean(x)
    return sum(x) / length(x)
end

function main()
    println("="^60)
    println("    SLAM: Scattering Lagrange Asymptotic Matching")
    println("    GUI using Blink.jl")
    println("="^60)
    println()

    # Create window with custom icon
    icon_path = joinpath(@__DIR__, "icon.png")

    # Create window with options
    global win = Window(Dict(
        :title => "SLAM Nuclear Scattering Calculator",
        :icon => icon_path
    ))
    size(win, 1800, 1000)

    # Load HTML
    body!(win, HTML_CONTENT)

    # Center window
    sleep(0.5)

    try
        win_width = 1800
        win_height = 1000

        Blink.js(win, Blink.JSString("""
        (function() {
            try {
                const screenWidth = window.screen.width;
                const screenHeight = window.screen.height;
                const screenX = window.screen.availLeft || 0;
                const screenY = window.screen.availTop || 0;

                const winWidth = $win_width;
                const winHeight = $win_height;
                const x = screenX + Math.round((screenWidth - winWidth) / 2);
                const y = screenY + Math.round((screenHeight - winHeight) / 2);

                window.moveTo(x, y);
            } catch (e) {
                console.error('Window positioning error:', e.message);
            }
        })();
        """))

        println("Window centered on screen")
    catch e
        println("Note: Could not apply window positioning: ", e)
    end

    # Load Plotly.js
    println("Loading Plotly.js...")
    plotly_path = joinpath(@__DIR__, "plotly-2.33.0.min.js")
    if isfile(plotly_path)
        plotly_js = read(plotly_path, String)
        Blink.js(win, Blink.JSString(plotly_js))
        println("Plotly.js loaded from local file")
    else
        # Load from CDN
        Blink.js(win, Blink.JSString("""
        var script = document.createElement('script');
        script.src = 'https://cdn.plot.ly/plotly-2.33.0.min.js';
        document.head.appendChild(script);
        """))
        println("Plotly.js loading from CDN...")
    end
    sleep(0.5)

    # Register message handler
    handle(win, "run_calculation") do params
        @async begin
            try
                run_calculation(params)
            catch e
                error_msg = sprint(showerror, e, catch_backtrace())
                println("Error in calculation: ", error_msg)
                try
                    full_error = "<br><br>ERROR:<br>" * error_msg
                    @js_ win document.getElementById("status").textContent = "Calculation error"
                    @js_ win document.getElementById("status").className = "status error"
                    @js_ win document.getElementById("output").innerHTML += $full_error
                    @js_ win document.getElementById("run-btn").disabled = false
                catch js_err
                    println("Could not send error to UI: ", js_err)
                end
            end
        end
        return nothing
    end

    println("\nGUI launched successfully!")
    println("  Adjust parameters and click 'Run All Partial Waves' to start.")
    println("  Close the window to exit.")
    println()

    # Keep window open
    try
        while active(win)
            sleep(0.1)
        end
    catch e
        if isa(e, InterruptException)
            println("\n\nShutting down...")
        else
            rethrow(e)
        end
    end
end

# Run the GUI
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
