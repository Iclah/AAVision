"""
HTML report generation for AAV analysis results.
Uses Jinja2 templating to create interactive reports with Plotly visualizations.
"""

import os
import json
import logging
from jinja2 import Template
from typing import Dict, List


def generate_html_report(stats: Dict, truncation_data: Dict, 
                        length_distributions: Dict, output_file: str,
                        template_path: str = None, coverage_data: Dict = None,
                        taxonomy_data: Dict = None) -> str:
    """
    Generate interactive HTML report from analysis results.
    
    Args:
        stats: Dictionary with classification statistics
        truncation_data: Dictionary mapping positions to enrichment percentages (rAAV only)
        length_distributions: Dict mapping sample names to length bins
        output_file: Path for output HTML file
        template_path: Path to Jinja2 template (optional)
        coverage_data: Dict mapping category names to coverage dicts
                      e.g., {'vector': {pos: depth}, 'helper': {...}, ...}
        taxonomy_data: Dict from taxonomy.run_taxonomy_analysis (optional)
                      Contains kraken_stats, kraken_report, bracken_results,
                      bracken_level, domain_summary
        
    Returns:
        Path to generated HTML report
    """
    logging.info("Generating HTML report")
    
    # Load template
    if template_path and os.path.exists(template_path):
        with open(template_path, 'r') as f:
            template_str = f.read()
    else:
        # Use embedded template
        template_str = get_default_template()
    
    template = Template(template_str)
    
    # Format numbers for display
    formatted_stats = format_statistics(stats)
    
    # Handle coverage_data
    if coverage_data is None:
        coverage_data = {}
    
    # Prepare data for JavaScript
    js_data = {
        # Truncation data (rAAV only)
        'truncation': json.dumps(truncation_data),
        # Coverage data for each category
        'coverage_vector': json.dumps(coverage_data.get('vector', {})),
        'coverage_helper': json.dumps(coverage_data.get('helper', {})) if 'helper' in coverage_data else None,
        'coverage_rep_cap': json.dumps(coverage_data.get('rep_cap', {})) if 'rep_cap' in coverage_data else None,
        'coverage_spike_in': json.dumps(coverage_data.get('spike_in', {})) if 'spike_in' in coverage_data else None,
        # Length distribution data
        'vector_lengths': json.dumps(length_distributions.get('vector', {})),
        'helper_lengths': json.dumps(length_distributions.get('helper', {})) if 'helper' in length_distributions else None,
        'rep_cap_lengths': json.dumps(length_distributions.get('rep_cap', {})) if 'rep_cap' in length_distributions else None,
        'host_lengths': json.dumps(length_distributions.get('host', {})) if 'host' in length_distributions else None,
        'spike_in_lengths': json.dumps(length_distributions.get('spike_in', {})) if 'spike_in' in length_distributions else None,
        'unmapped_lengths': json.dumps(length_distributions.get('unmapped', {})) if 'unmapped' in length_distributions else None,
        # Taxonomy data
        'taxonomy_bracken': json.dumps(taxonomy_data.get('bracken_results', [])) if taxonomy_data else 'null',
        'taxonomy_domains': json.dumps(taxonomy_data.get('domain_summary', {})) if taxonomy_data else 'null',
    }
    
    # Determine which optional tabs to show for read length (includes host)
    optional_tabs = []
    if 'helper' in length_distributions and length_distributions['helper']:
        optional_tabs.append('helper')
    if 'rep_cap' in length_distributions and length_distributions['rep_cap']:
        optional_tabs.append('rep_cap')
    if 'host' in length_distributions and length_distributions['host']:
        optional_tabs.append('host')
    if 'spike_in' in length_distributions and length_distributions['spike_in']:
        optional_tabs.append('spike_in')
    if 'unmapped' in length_distributions and length_distributions['unmapped']:
        optional_tabs.append('unmapped')
    
    # Determine which coverage tabs to show (excludes host)
    coverage_tabs = ['vector']  # Always have vector
    if 'helper' in coverage_data and coverage_data['helper']:
        coverage_tabs.append('helper')
    if 'rep_cap' in coverage_data and coverage_data['rep_cap']:
        coverage_tabs.append('rep_cap')
    if 'spike_in' in coverage_data and coverage_data['spike_in']:
        coverage_tabs.append('spike_in')
    
    # Render template
    html_content = template.render(
        stats=formatted_stats,
        js_data=js_data,
        optional_tabs=optional_tabs,
        coverage_tabs=coverage_tabs,
        taxonomy_data=taxonomy_data
    )
    
    # Write to file
    with open(output_file, 'w') as f:
        f.write(html_content)
    
    logging.info(f"HTML report generated: {output_file}")
    return output_file


def format_statistics(stats: Dict) -> Dict:
    """
    Format statistics for display in HTML report.
    
    Args:
        stats: Raw statistics dictionary
        
    Returns:
        Formatted statistics dictionary
    """
    formatted = stats.copy()
    
    # Format large numbers with commas
    for key in ['total_reads', 'total_vector', 'total_hybrid', 'total_impurity',
                'vector_hq', 'vector_lq_pure', 'backbone_hybrid', 'helper_hybrid',
                'host_hybrid', 'backbone_impurity', 'helper_impurity', 'host_impurity']:
        if key in formatted:
            formatted[f'{key}_formatted'] = f"{formatted[key]:,}"
    
    # Format percentages
    for key in ['vector_percent', 'hybrid_percent', 'impurity_percent']:
        if key in formatted:
            formatted[f'{key}_rounded'] = f"{formatted[key]:.2f}"
    
    # Calculate individual percentages
    total = formatted.get('total_reads', 1)
    if total > 0:
        for key in ['vector_hq', 'vector_lq_pure', 'backbone_hybrid', 'helper_hybrid',
                    'host_hybrid', 'backbone_impurity', 'helper_impurity', 'host_impurity']:
            if key in formatted:
                formatted[f'{key}_percent'] = f"{(formatted[key] / total * 100):.2f}"
    
    # Handle optional references
    for ref in ['rep_cap', 'spike_in']:
        for category in ['hybrid', 'impurity']:
            key = f'{ref}_{category}'
            if key in formatted:
                formatted[f'{key}_formatted'] = f"{formatted[key]:,}"
                if total > 0:
                    formatted[f'{key}_percent'] = f"{(formatted[key] / total * 100):.2f}"
    
    return formatted


def get_default_template() -> str:
    """
    Return the default HTML template as a string.
    This is embedded to avoid external file dependencies.
    
    Returns:
        HTML template string
    """
    return '''<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>AAVision Analysis Report</title>
  <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
  <link href="https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700&display=swap" rel="stylesheet">
  <style>
    * {
      margin: 0;
      padding: 0;
      box-sizing: border-box;
    }
    
    body {
      background: linear-gradient(135deg, #1e1e2f 0%, #16213e 100%);
      color: #e8e8e8;
      font-family: 'Inter', sans-serif;
      min-height: 100vh;
      padding: 20px;
      line-height: 1.6;
    }
    
    .header {
      text-align: center;
      margin-bottom: 40px;
    }
    
    h1 {
      font-size: 2.8rem;
      font-weight: 700;
      background: linear-gradient(135deg, #4facfe 0%, #00f2fe 100%);
      -webkit-background-clip: text;
      -webkit-text-fill-color: transparent;
      background-clip: text;
      margin-bottom: 10px;
    }
    
    .subtitle {
      color: #a0a0a0;
      font-size: 1.1rem;
      font-weight: 300;
    }
    
    .stats-overview {
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
      gap: 20px;
      margin-bottom: 30px;
      max-width: 1200px;
      margin-left: auto;
      margin-right: auto;
    }
    
    .stat-card {
      background: rgba(255, 255, 255, 0.05);
      backdrop-filter: blur(10px);
      border: 1px solid rgba(255, 255, 255, 0.1);
      border-radius: 16px;
      padding: 20px;
      text-align: center;
      transition: all 0.3s ease;
    }
    
    .stat-card:hover {
      transform: translateY(-5px);
      box-shadow: 0 20px 40px rgba(0, 0, 0, 0.3);
      border-color: rgba(255, 255, 255, 0.2);
    }
    
    .stat-value {
      font-size: 2.5rem;
      font-weight: 700;
      margin-bottom: 5px;
    }
    
    .stat-label {
      color: #a0a0a0;
      font-size: 0.9rem;
      letter-spacing: 1px;
    }
    
    .raav { color: #4caf50; }
    .hybrid { color: #ff9800; }
    .impurity { color: #f44336; }
    
    .card {
      background: rgba(255, 255, 255, 0.05);
      backdrop-filter: blur(20px);
      border: 1px solid rgba(255, 255, 255, 0.1);
      border-radius: 20px;
      padding: 30px;
      margin: 30px auto;
      max-width: 1200px;
      box-shadow: 0 15px 35px rgba(0, 0, 0, 0.2);
      transition: all 0.3s ease;
    }
    
    .card:hover {
      border-color: rgba(255, 255, 255, 0.2);
    }
    
    .card h2 {
      font-size: 1.8rem;
      font-weight: 600;
      margin-bottom: 25px;
      color: #f0f0f0;
    }
    
    .distribution-container {
      display: grid;
      grid-template-columns: 1fr 1fr;
      gap: 40px;
      align-items: start;
    }
    
    table {
      width: 100%;
      border-collapse: collapse;
    }
    
    th, td {
      padding: 12px 16px;
      text-align: right;
      border-bottom: 1px solid rgba(255, 255, 255, 0.1);
      transition: all 0.2s ease;
    }
    
    th:first-child, td:first-child {
      text-align: left;
      width: 60%;
    }
    
    th {
      background: rgba(255, 255, 255, 0.1);
      font-weight: 600;
      color: #f0f0f0;
      text-transform: uppercase;
      font-size: 0.85rem;
      letter-spacing: 1px;
    }
    
    tr:hover {
      background: rgba(255, 255, 255, 0.05);
    }
    
    .category-row {
      font-weight: 600;
      font-size: 1.05rem;
    }
    
    .raav-row { background: rgba(76, 175, 80, 0.4); }
    .hybrid-row { background: rgba(255, 152, 0, 0.4); }
    .impurity-row { background: rgba(244, 67, 54, 0.4); }
    
    .subcategory td {
      padding-left: 32px;
      color: #b0b0b0;
      font-size: 0.95rem;
    }
    
    .chart-container {
      height: 440px;
      position: relative;
      padding: 40px 20px 20px 20px;
    }
    
    .tabs {
      display: flex;
      border-bottom: 2px solid rgba(255, 255, 255, 0.1);
      margin-bottom: 25px;
    }
    
    .tab {
      padding: 12px 24px;
      background: transparent;
      border: none;
      color: #a0a0a0;
      font-family: 'Inter', sans-serif;
      font-size: 0.95rem;
      font-weight: 500;
      cursor: pointer;
      transition: all 0.3s ease;
      border-radius: 8px 8px 0 0;
      margin-right: 4px;
    }
    
    .tab:hover {
      color: #e0e0e0;
      background: rgba(255, 255, 255, 0.05);
    }
    
    .tab.active {
      color: #4facfe;
      background: rgba(79, 172, 254, 0.1);
      border-bottom: 2px solid #4facfe;
    }
    
    .tab-content {
      display: none;
      height: 450px;
    }
    
    .tab-content.active {
      display: block;
    }
    
    @media (max-width: 768px) {
      .distribution-container {
        grid-template-columns: 1fr;
        gap: 30px;
      }
      
      h1 {
        font-size: 2.2rem;
      }
      
      .stats-overview {
        grid-template-columns: 1fr;
      }
      
      .card {
        padding: 20px;
      }
    }
  </style>
</head>
<body>
  <div class="header">
    <h1>AAVision Analysis Report</h1>
    <p class="subtitle">Comprehensive AAV sequencing analysis</p>
  </div>

  <div class="stats-overview">
    <div class="stat-card">
      <div class="stat-value raav">{{ stats.vector_percent_rounded }}%</div>
      <div class="stat-label">rAAV</div>
    </div>
    <div class="stat-card">
      <div class="stat-value hybrid">{{ stats.hybrid_percent_rounded }}%</div>
      <div class="stat-label">Hybrid</div>
    </div>
    <div class="stat-card">
      <div class="stat-value impurity">{{ stats.impurity_percent_rounded }}%</div>
      <div class="stat-label">Impurity</div>
    </div>
    <div class="stat-card">
      <div class="stat-value">{{ stats.total_reads_formatted }}</div>
      <div class="stat-label">Total Reads</div>
    </div>
  </div>

  <div class="card">
    <h2>Read Distribution</h2>
    <div class="distribution-container">
      <div class="table-container">
        <table>
          <thead>
            <tr><th>Category</th><th>Reads</th><th>Percentage</th></tr>
          </thead>
          <tbody>
            <tr class="category-row raav-row">
              <td><strong>rAAV</strong></td>
              <td>{{ stats.total_vector_formatted }}</td>
              <td><strong>{{ stats.vector_percent_rounded }}%</strong></td>
            </tr>
            <tr class="subcategory">
              <td>High Quality</td>
              <td>{{ stats.vector_hq_formatted }}</td>
              <td>{{ stats.vector_hq_percent }}%</td>
            </tr>
            <tr class="subcategory">
              <td>Low Quality</td>
              <td>{{ stats.vector_lq_pure_formatted }}</td>
              <td>{{ stats.vector_lq_pure_percent }}%</td>
            </tr>
            <tr class="category-row hybrid-row">
              <td><strong>Hybrid</strong></td>
              <td>{{ stats.total_hybrid_formatted }}</td>
              <td><strong>{{ stats.hybrid_percent_rounded }}%</strong></td>
            </tr>
            <tr class="subcategory">
              <td>Backbone Hybrid</td>
              <td>{{ stats.backbone_hybrid_formatted }}</td>
              <td>{{ stats.backbone_hybrid_percent }}%</td>
            </tr>
            <tr class="subcategory">
              <td>Helper Hybrid</td>
              <td>{{ stats.helper_hybrid_formatted }}</td>
              <td>{{ stats.helper_hybrid_percent }}%</td>
            </tr>
            <tr class="subcategory">
              <td>Host Genome Hybrid</td>
              <td>{{ stats.host_hybrid_formatted }}</td>
              <td>{{ stats.host_hybrid_percent }}%</td>
            </tr>
            {% if stats.rep_cap_hybrid_formatted %}
            <tr class="subcategory">
              <td>Rep/Cap Hybrid</td>
              <td>{{ stats.rep_cap_hybrid_formatted }}</td>
              <td>{{ stats.rep_cap_hybrid_percent }}%</td>
            </tr>
            {% endif %}
            <tr class="category-row impurity-row">
              <td><strong>Impurity</strong></td>
              <td>{{ stats.total_impurity_formatted }}</td>
              <td><strong>{{ stats.impurity_percent_rounded }}%</strong></td>
            </tr>
            <tr class="subcategory">
              <td>Backbone Impurity</td>
              <td>{{ stats.backbone_impurity_formatted }}</td>
              <td>{{ stats.backbone_impurity_percent }}%</td>
            </tr>
            <tr class="subcategory">
              <td>Helper Impurity</td>
              <td>{{ stats.helper_impurity_formatted }}</td>
              <td>{{ stats.helper_impurity_percent }}%</td>
            </tr>
            <tr class="subcategory">
              <td>Host Genome Impurity</td>
              <td>{{ stats.host_impurity_formatted }}</td>
              <td>{{ stats.host_impurity_percent }}%</td>
            </tr>
            {% if stats.rep_cap_impurity_formatted %}
            <tr class="subcategory">
              <td>Rep/Cap Impurity</td>
              <td>{{ stats.rep_cap_impurity_formatted }}</td>
              <td>{{ stats.rep_cap_impurity_percent }}%</td>
            </tr>
            {% endif %}
            {% if stats.spike_in_impurity_formatted %}
            <tr class="subcategory">
              <td>Spike-in Impurity</td>
              <td>{{ stats.spike_in_impurity_formatted }}</td>
              <td>{{ stats.spike_in_impurity_percent }}%</td>
            </tr>
            {% endif %}
          </tbody>
        </table>
      </div>
      <div class="chart-container">
        <div id="pieChart"></div>
      </div>
    </div>
  </div>

  <div class="card">
    <h2>Read Length Distribution</h2>
    <div class="tabs">
      <button class="tab active" onclick="showTab(event, 'raav-hist')">rAAV Reads</button>
      {% if 'helper' in optional_tabs %}
      <button class="tab" onclick="showTab(event, 'helper-hist')">Helper Impurity</button>
      {% endif %}
      {% if 'rep_cap' in optional_tabs %}
      <button class="tab" onclick="showTab(event, 'repcap-hist')">Rep/Cap Impurity</button>
      {% endif %}
      {% if 'host' in optional_tabs %}
      <button class="tab" onclick="showTab(event, 'host-hist')">Host Genome Impurity</button>
      {% endif %}
      {% if 'spike_in' in optional_tabs %}
      <button class="tab" onclick="showTab(event, 'spikein-hist')">Spike-in Reads</button>
      {% endif %}
      {% if 'unmapped' in optional_tabs %}
      <button class="tab" onclick="showTab(event, 'unmapped-hist')">Unclassified</button>
      {% endif %}
    </div>
    
    <div id="raav-hist" class="tab-content active">
      <div id="histPlot"></div>
    </div>
    {% if 'helper' in optional_tabs %}
    <div id="helper-hist" class="tab-content">
      <div id="helperHistPlot"></div>
    </div>
    {% endif %}
    {% if 'rep_cap' in optional_tabs %}
    <div id="repcap-hist" class="tab-content">
      <div id="repCapHistPlot"></div>
    </div>
    {% endif %}
    {% if 'host' in optional_tabs %}
    <div id="host-hist" class="tab-content">
      <div id="hostHistPlot"></div>
    </div>
    {% endif %}
    {% if 'spike_in' in optional_tabs %}
    <div id="spikein-hist" class="tab-content">
      <div id="spikeInHistPlot"></div>
    </div>
    {% endif %}
    {% if 'unmapped' in optional_tabs %}
    <div id="unmapped-hist" class="tab-content">
      <div id="unmappedHistPlot"></div>
    </div>
    {% endif %}
  </div>

  <div class="card">
    <h2>Truncation Hotspots</h2>
    <div class="chart-container">
      <div id="truncationPlot"></div>
    </div>
  </div>

  <div class="card">
    <h2>Read Coverage</h2>
    <div class="tabs">
      <button class="tab active" onclick="showTab(event, 'cov-raav')">rAAV Reads</button>
      {% if 'helper' in coverage_tabs %}
      <button class="tab" onclick="showTab(event, 'cov-helper')">Helper Impurity</button>
      {% endif %}
      {% if 'rep_cap' in coverage_tabs %}
      <button class="tab" onclick="showTab(event, 'cov-repcap')">Rep/Cap Impurity</button>
      {% endif %}
      {% if 'spike_in' in coverage_tabs %}
      <button class="tab" onclick="showTab(event, 'cov-spikein')">Spike-in Reads</button>
      {% endif %}
    </div>
    
    <div id="cov-raav" class="tab-content active">
      <div id="coveragePlot"></div>
    </div>
    {% if 'helper' in coverage_tabs %}
    <div id="cov-helper" class="tab-content">
      <div id="coverageHelperPlot"></div>
    </div>
    {% endif %}
    {% if 'rep_cap' in coverage_tabs %}
    <div id="cov-repcap" class="tab-content">
      <div id="coverageRepCapPlot"></div>
    </div>
    {% endif %}
    {% if 'spike_in' in coverage_tabs %}
    <div id="cov-spikein" class="tab-content">
      <div id="coverageSpikeInPlot"></div>
    </div>
    {% endif %}
  </div>

  {% if taxonomy_data %}
  <div class="card">
    <h2>Taxonomic Classification of Unmapped Reads</h2>
    <p style="color: #a0a0a0; margin-bottom: 20px; font-size: 0.9rem;">
      Kraken2 classification of reads that did not map to any provided reference.
      {% if taxonomy_data.bracken_results %}Abundance re-estimated with Bracken at {{ taxonomy_data.bracken_level }}-level.{% endif %}
    </p>
    
    <div class="distribution-container">
      <div class="table-container">
        {% if taxonomy_data.bracken_results %}
        <table>
          <thead>
            <tr><th>Taxon</th><th>Reads</th><th>Abundance</th></tr>
          </thead>
          <tbody>
            {% for entry in taxonomy_data.bracken_results %}
            <tr>
              <td style="font-style: italic;">{{ entry.name }}</td>
              <td>{{ "{:,}".format(entry.reads) }}</td>
              <td>{{ entry.fraction_percent }}%</td>
            </tr>
            {% endfor %}
          </tbody>
        </table>
        {% else %}
        <table>
          <thead>
            <tr><th>Taxon</th><th>Reads (clade)</th><th>Percentage</th></tr>
          </thead>
          <tbody>
            {% for entry in taxonomy_data.kraken_report %}
            <tr>
              <td>
                <span style="color: #666; font-size: 0.8rem;">[{{ entry.rank }}]</span>
                {% if entry.rank in ['S', 'G'] %}<span style="font-style: italic;">{{ entry.name }}</span>{% else %}{{ entry.name }}{% endif %}
              </td>
              <td>{{ "{:,}".format(entry.reads_clade) }}</td>
              <td>{{ entry.percent }}%</td>
            </tr>
            {% endfor %}
          </tbody>
        </table>
        {% endif %}
      </div>
      <div class="chart-container">
        <div id="taxonomyPieChart"></div>
      </div>
    </div>
    
    <div style="margin-top: 20px; padding: 15px; background: rgba(255,255,255,0.03); border-radius: 10px; font-size: 0.85rem; color: #888;">
      Kraken2: {{ "{:,}".format(taxonomy_data.kraken_stats.classified) }} classified ({{ taxonomy_data.kraken_stats.classified_percent }}%),
      {{ "{:,}".format(taxonomy_data.kraken_stats.unclassified) }} unclassified ({{ taxonomy_data.kraken_stats.unclassified_percent }}%)
      of {{ "{:,}".format(taxonomy_data.kraken_stats.total) }} input reads.
    </div>
  </div>
  {% endif %}

  <script>
    const plotConfig = {
      responsive: true,
      displayModeBar: true,
      modeBarButtonsToRemove: ['pan2d', 'lasso2d', 'select2d', 'autoScale2d', 'resetScale2d'],
      toImageButtonOptions: {
        format: 'png',
        filename: 'aav_analysis_chart',
        height: 600,
        width: 800,
        scale: 2
      }
    };

    const darkLayout = {
      paper_bgcolor: 'rgba(0,0,0,0)',
      plot_bgcolor: 'rgba(0,0,0,0)',
      font: { color: '#ccc', family: 'Inter, sans-serif' },
      hoverlabel: {
        bgcolor: 'rgba(50, 50, 60, 0.9)',
        bordercolor: 'rgba(255,255,255,0.2)',
        font: { color: '#fff' }
      }
    };

    // Helper function to create histogram
    function createHistogram(elementId, bins, color, hoverColor) {
      if (!bins || Object.keys(bins).length === 0) return;
      
      const labels = Object.keys(bins).sort((a, b) => {
        // Handle overflow bin (≥6000 or >6000) - sort to end
        const aIsOverflow = a.startsWith('≥') || a.startsWith('>');
        const bIsOverflow = b.startsWith('≥') || b.startsWith('>');
        if (aIsOverflow && !bIsOverflow) return 1;
        if (!aIsOverflow && bIsOverflow) return -1;
        if (aIsOverflow && bIsOverflow) return 0;
        // Normal bins - sort by start number
        const aNum = parseInt(a.split('-')[0]);
        const bNum = parseInt(b.split('-')[0]);
        return aNum - bNum;
      });
      const percentages = labels.map(label => bins[label]);
      const maxPercentage = Math.max(...percentages);
      const yAxisMax = Math.ceil(maxPercentage * 1.1) || 1;
      
      const histData = [{
        x: labels,
        y: percentages,
        type: 'bar',
        marker: { 
          color: color,
          line: { color: 'rgba(255,255,255,0.1)', width: 1 }
        },
        hovertemplate: '<b>%{x} bp</b><br>Percentage: %{y:.2f}%<br><extra></extra>',
        hoverlabel: { bgcolor: hoverColor }
      }];

      const histLayout = {
        ...darkLayout,
        yaxis: { 
          title: { text: 'Percentage (%)', font: { color: '#ccc' } },
          gridcolor: '#444',
          zerolinecolor: '#666',
          range: [0, yAxisMax]
        },
        xaxis: { 
          title: { text: 'Read Length (bp)', font: { color: '#ccc' } },
          tickangle: -45,
          zerolinecolor: '#666'
        },
        margin: { t: 20, b: 80, l: 60, r: 20 }
      };

      Plotly.newPlot(elementId, histData, histLayout, plotConfig);
    }

    // Helper function to create coverage plot
    function createCoveragePlot(elementId, data, color, hoverColor) {
      if (!data || Object.keys(data).length === 0) {
        document.getElementById(elementId).innerHTML = '<p style="text-align: center; color: #a0a0a0; padding: 50px;">No coverage data available</p>';
        return;
      }
      
      const positions = Object.keys(data).map(pos => parseInt(pos)).sort((a, b) => a - b);
      const depths = positions.map(pos => data[pos.toString()]);
      const maxDepth = Math.max(...depths);
      const yAxisMax = Math.ceil(maxDepth * 1.1) || 1;

      const coveragePlotData = [{
        x: positions,
        y: depths,
        type: 'scatter',
        mode: 'lines',
        name: 'Coverage Depth',
        line: { color: color, width: 2 },
        fill: 'tozeroy',
        fillcolor: color.replace('rgb', 'rgba').replace(')', ', 0.3)'),
        hovertemplate: '<b>Position: %{x} bp</b><br>Coverage: %{y:,}<br><extra></extra>',
        hoverlabel: { bgcolor: hoverColor }
      }];

      const coverageLayout = {
        ...darkLayout,
        yaxis: {
          title: { text: 'Coverage Depth (reads)', font: { color: '#ccc', size: 14 } },
          gridcolor: 'rgba(255,255,255,0.1)',
          zerolinecolor: 'rgba(255,255,255,0.2)',
          range: [0, yAxisMax]
        },
        xaxis: {
          title: { text: 'Genomic Position (bp)', font: { color: '#ccc', size: 14 } },
          gridcolor: 'rgba(255,255,255,0.1)',
          zerolinecolor: 'rgba(255,255,255,0.2)',
          dtick: 500,
          tickmode: 'linear',
          tickangle: -45
        },
        margin: { t: 60, b: 80, l: 80, r: 40 },
        showlegend: false
      };

      Plotly.newPlot(elementId, coveragePlotData, coverageLayout, plotConfig);
    }

    // Pie Chart
    const pieData = [{
      values: [{{ stats.total_vector }}, {{ stats.total_hybrid }}, {{ stats.total_impurity }}],
      labels: ['rAAV', 'Hybrid', 'Impurity'],
      type: 'pie',
      marker: { 
        colors: ['#4caf50', '#ff9800', '#f44336'],
        line: { color: 'rgba(255,255,255,0.1)', width: 2 }
      },
      textinfo: 'label+percent',
      textposition: 'auto',
      textfont: { color: '#ffffff', size: 14, family: 'Inter' },
      hovertemplate: '<b>%{label}</b><br>Reads: %{value:,}<br>Percentage: %{percent}<br><extra></extra>',
      pull: [0, 0, 0]
    }];

    const pieLayout = {
      ...darkLayout,
      margin: { t: 80, b: 50, l: 50, r: 50 },
      showlegend: false
    };

    Plotly.newPlot('pieChart', pieData, pieLayout, plotConfig);

    // Create histograms
    const vectorBins = {{ js_data.vector_lengths|safe }};
    createHistogram('histPlot', vectorBins, '#4caf50', 'rgba(76, 175, 80, 0.8)');

    {% if js_data.helper_lengths %}
    const helperBins = {{ js_data.helper_lengths|safe }};
    createHistogram('helperHistPlot', helperBins, '#9c27b0', 'rgba(156, 39, 176, 0.8)');
    {% endif %}

    {% if js_data.rep_cap_lengths %}
    const repCapBins = {{ js_data.rep_cap_lengths|safe }};
    createHistogram('repCapHistPlot', repCapBins, '#2196f3', 'rgba(33, 150, 243, 0.8)');
    {% endif %}

    {% if js_data.host_lengths %}
    const hostBins = {{ js_data.host_lengths|safe }};
    createHistogram('hostHistPlot', hostBins, '#ff5722', 'rgba(255, 87, 34, 0.8)');
    {% endif %}

    {% if js_data.spike_in_lengths %}
    const spikeInBins = {{ js_data.spike_in_lengths|safe }};
    createHistogram('spikeInHistPlot', spikeInBins, '#f44336', 'rgba(244, 67, 54, 0.8)');
    {% endif %}

    {% if js_data.unmapped_lengths %}
    const unmappedBins = {{ js_data.unmapped_lengths|safe }};
    createHistogram('unmappedHistPlot', unmappedBins, '#9e9e9e', 'rgba(158, 158, 158, 0.8)');
    {% endif %}

    // Truncation plot (rAAV only)
    const truncationData = {{ js_data.truncation|safe }};
    
    if (Object.keys(truncationData).length > 0) {
      const positions = Object.keys(truncationData).map(pos => parseInt(pos)).sort((a, b) => a - b);
      const percentages = positions.map(pos => truncationData[pos.toString()]);
      const maxPercentage = Math.max(...percentages);
      const yAxisMax = Math.ceil(maxPercentage * 1.1) || 1;
      
      const truncationPlotData = [{
        x: positions,
        y: percentages,
        type: 'scatter',
        mode: 'lines',
        name: 'Enrichment Score',
        line: { color: '#00f2fe', width: 4 },
        hovertemplate: '<b>Position: %{x} bp</b><br>Enrichment: %{y:.2f}%<br><extra></extra>',
        hoverlabel: { bgcolor: 'rgba(0, 242, 254, 0.8)' }
      }];

      const truncationLayout = {
        ...darkLayout,
        yaxis: { 
          title: { text: 'Enrichment Score (%)', font: { color: '#ccc', size: 14 } },
          gridcolor: 'rgba(255,255,255,0.1)',
          zerolinecolor: 'rgba(255,255,255,0.2)',
          range: [0, yAxisMax]
        },
        xaxis: { 
          title: { text: 'Genomic Position (bp)', font: { color: '#ccc', size: 14 } },
          gridcolor: 'rgba(255,255,255,0.1)',
          zerolinecolor: 'rgba(255,255,255,0.2)',
          dtick: 500,
          tickmode: 'linear',
          tickangle: -45
        },
        margin: { t: 60, b: 80, l: 80, r: 40 },
        showlegend: false
      };

      Plotly.newPlot('truncationPlot', truncationPlotData, truncationLayout, plotConfig);
    }

    // Create coverage plots
    const coverageVectorData = {{ js_data.coverage_vector|safe }};
    createCoveragePlot('coveragePlot', coverageVectorData, 'rgb(76, 175, 80)', 'rgba(76, 175, 80, 0.8)');

    {% if js_data.coverage_helper %}
    const coverageHelperData = {{ js_data.coverage_helper|safe }};
    createCoveragePlot('coverageHelperPlot', coverageHelperData, 'rgb(156, 39, 176)', 'rgba(156, 39, 176, 0.8)');
    {% endif %}

    {% if js_data.coverage_rep_cap %}
    const coverageRepCapData = {{ js_data.coverage_rep_cap|safe }};
    createCoveragePlot('coverageRepCapPlot', coverageRepCapData, 'rgb(33, 150, 243)', 'rgba(33, 150, 243, 0.8)');
    {% endif %}

    {% if js_data.coverage_spike_in %}
    const coverageSpikeInData = {{ js_data.coverage_spike_in|safe }};
    createCoveragePlot('coverageSpikeInPlot', coverageSpikeInData, 'rgb(244, 67, 54)', 'rgba(244, 67, 54, 0.8)');
    {% endif %}

    function showTab(event, tabId) {
      // Find parent card to scope tab switching
      const card = event.currentTarget.closest('.card');
      card.querySelectorAll('.tab-content').forEach(content => content.classList.remove('active'));
      card.querySelectorAll('.tab').forEach(tab => tab.classList.remove('active'));
      card.querySelector('#' + tabId).classList.add('active');
      event.currentTarget.classList.add('active');
      
      setTimeout(() => {
        const activeTabContent = card.querySelector('.tab-content.active');
        const plotDiv = activeTabContent.querySelector('[id*="Plot"]');
        if (plotDiv) {
          Plotly.Plots.resize(plotDiv);
        }
      }, 100);
    }

    // Pie chart hover animations
    const pieChartElement = document.getElementById('pieChart');
    
    pieChartElement.on('plotly_hover', function(eventData) {
      const pointIndex = eventData.points[0].pointNumber;
      const currentPull = [0, 0, 0];
      currentPull[pointIndex] = 0.15;
      Plotly.restyle('pieChart', { 'pull': [currentPull] });
    });

    pieChartElement.on('plotly_unhover', function(eventData) {
      Plotly.restyle('pieChart', { 'pull': [[0, 0, 0]] });
    });

    // --- Taxonomy pie chart ---
    {% if taxonomy_data %}
    const taxonomyDomains = {{ js_data.taxonomy_domains|safe }};
    const taxonomyBracken = {{ js_data.taxonomy_bracken|safe }};
    
    // Use Bracken results if available, otherwise use domain summary
    let taxLabels, taxValues, taxColors;
    
    if (taxonomyBracken && taxonomyBracken.length > 0) {
      // Show top Bracken species/genus results
      const topN = taxonomyBracken.slice(0, 10);
      const topReads = topN.reduce((sum, e) => sum + e.reads, 0);
      const topSumPercent = topN.reduce((sum, e) => sum + e.fraction_percent, 0);
      
      taxLabels = topN.map(e => e.name);
      taxValues = topN.map(e => e.reads);
      
      // If top N doesn't cover ~100%, add an "Other" slice using implied read count
      if (topSumPercent < 99.5 && topSumPercent > 0) {
        const impliedTotal = topReads / (topSumPercent / 100);
        const otherReads = Math.max(0, Math.round(impliedTotal - topReads));
        if (otherReads > 0) {
          taxLabels.push('Other');
          taxValues.push(otherReads);
        }
      }
      
      // Generate colors from a palette
      const palette = [
        '#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231',
        '#911eb4', '#42d4f4', '#f032e6', '#bfef45', '#fabed4', '#9e9e9e'
      ];
      taxColors = taxLabels.map((_, i) => palette[i % palette.length]);
    } else if (taxonomyDomains && Object.keys(taxonomyDomains).length > 0) {
      // Domain-level summary fallback
      taxLabels = Object.keys(taxonomyDomains);
      taxValues = taxLabels.map(d => taxonomyDomains[d].reads);
      
      const domainColors = {
        'Bacteria': '#e6194b',
        'Viruses': '#3cb44b',
        'Eukaryota': '#4363d8',
        'Archaea': '#f58231',
        'Unclassified': '#9e9e9e'
      };
      taxColors = taxLabels.map(d => domainColors[d] || '#888888');
    } else {
      taxLabels = [];
      taxValues = [];
      taxColors = [];
    }
    
    if (taxLabels.length > 0) {
      const taxPieData = [{
        values: taxValues,
        labels: taxLabels,
        type: 'pie',
        marker: {
          colors: taxColors,
          line: { color: 'rgba(255,255,255,0.1)', width: 2 }
        },
        textinfo: 'label+percent',
        textposition: 'auto',
        textfont: { color: '#ffffff', size: 12, family: 'Inter' },
        hovertemplate: '<b>%{label}</b><br>Reads: %{value:,}<br>Percentage: %{percent}<br><extra></extra>'
      }];
      
      const taxPieLayout = {
        ...darkLayout,
        margin: { t: 40, b: 40, l: 40, r: 40 },
        showlegend: false
      };
      
      Plotly.newPlot('taxonomyPieChart', taxPieData, taxPieLayout, plotConfig);
    }
    {% endif %}
  </script>
</body>
</html>'''
