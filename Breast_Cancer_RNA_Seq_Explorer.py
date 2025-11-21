import os
import sys

# Fix for Dash Jupyter comm error when running outside Jupyter
# Patch the comm module before importing dash to avoid NotImplementedError
if 'comm' not in sys.modules:
    sys.modules['comm'] = type(sys)('comm')
    sys.modules['comm'].create_comm = lambda *a, **k: None

import dash
from dash import dcc, html, Input, Output, State
import plotly.graph_objs as go
import pandas as pd
import base64
import io

# Add the backend directory to the path to import backend_api
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
BACKEND_DIR = os.path.join(BASE_DIR, "backend")
if BACKEND_DIR not in sys.path:
    sys.path.insert(0, BACKEND_DIR)

# Import backend API functions
try:
    from backend_api import get_survival_analysis, get_gene_importance, get_expression_by_subtype, TCGA_DF, SUBTYPE_MAP, IMPORTANCE_DF
    BACKEND_AVAILABLE = True
except Exception as e:
    print(f"Warning: Could not import backend API: {e}")
    BACKEND_AVAILABLE = False
    TCGA_DF = None
    SUBTYPE_MAP = {}
    IMPORTANCE_DF = pd.DataFrame()

# Initialize Dash app
app = dash.Dash(__name__)
app.title = "Breast Cancer RNA Seq Explorer"

# Load actual gene lists from datasets
def get_available_genes(dataset='TCGA'):
    """Get list of available genes from the dataset"""
    df = TCGA_DF
    if df is None:
        return ['TP53', 'BRCA1', 'BRCA2', 'ESR1', 'HER2']  # Fallback
    # Exclude non-gene columns
    exclude_cols = ['time', 'status', 'subtype', 'PATIENT_ID']
    genes = [col for col in df.columns if col not in exclude_cols]
    return sorted(genes)

def get_available_subtypes(dataset='TCGA'):
    """Get list of available subtypes for the dataset"""
    df = TCGA_DF
    if df is None:
        return ['Luminal A', 'Luminal B', 'Basal-like', 'HER2-enriched', 'Normal-like']  # Fallback
    
    # Get unique subtypes and map to human-readable names
    unique_subtypes = df['subtype'].dropna().unique().tolist()  # Filter out NaN values
    # Reverse map from internal codes to human-readable
    reverse_map = {}
    if dataset in SUBTYPE_MAP:
        for human_name, internal_code in SUBTYPE_MAP[dataset].items():
            if internal_code in unique_subtypes:
                reverse_map[internal_code] = human_name
    
    # Return human-readable names, or internal codes if no mapping found
    readable = []
    for subtype in unique_subtypes:
        # Skip NaN, None, or empty values
        if pd.isna(subtype) or subtype == '':
            continue
        # Check if this subtype has a human-readable mapping
        if subtype in reverse_map:
            readable.append(reverse_map[subtype])
        else:
            # Use string representation for consistency
            readable.append(str(subtype))
    
    # Filter out any remaining None/NaN and ensure all are strings, then sort
    readable = [str(s) for s in readable if s and not pd.isna(s) and str(s).strip()]
    return sorted(set(readable))

# Initialize with TCGA data - All genes available (searchable dropdown handles large lists efficiently)
GENES = get_available_genes('TCGA') if BACKEND_AVAILABLE else ['TP53', 'BRCA1', 'BRCA2', 'ESR1', 'HER2']
SUBTYPES = get_available_subtypes('TCGA') if BACKEND_AVAILABLE else ['Luminal A', 'Luminal B', 'Basal-like', 'HER2-enriched', 'Normal-like']

# Get list of top 2000 genes (those with SHAP scores)
def get_top_genes():
    """Returns list of genes that have SHAP scores (top 2000)"""
    if not BACKEND_AVAILABLE:
        return []
    try:
        if not IMPORTANCE_DF.empty:
            return IMPORTANCE_DF.index.tolist()
    except:
        pass
    return []

TOP_GENES = get_top_genes()

# Layout - Centered and streamlined
app.layout = html.Div([
    html.Div([
        html.H1("🧬 Breast Cancer Survival Explorer", 
                style={'textAlign': 'center', 'marginBottom': '30px', 'color': '#2c3e50'}),
        
        # Control Panel - Centered
        html.Div([
            html.Div([
                html.Label("Search Gene", style={'fontWeight': 'bold', 'marginBottom': '5px'}),
                dcc.Checklist(
                    id='filter-top-genes',
                    options=[{'label': ' Show only top 2000 genes (with SHAP scores)', 'value': 'filter'}],
                    value=[],
                    style={'marginBottom': '10px', 'fontSize': '14px'}
                ),
                dcc.Dropdown(
                    id='gene-selector',
                    options=[{'label': gene, 'value': gene} for gene in GENES],
                    multi=False,
                    value='TP53',
                    searchable=True,
                    placeholder="Type to search for a gene...",
                    style={'marginBottom': '20px'}
                ),

                html.Label("Select Subtype(s) - Optional (leave empty for all subtypes)", 
                          style={'fontWeight': 'bold', 'marginBottom': '5px'}),
                dcc.Checklist(
                    id='subtype-selector',
                    options=[{'label': subtype, 'value': subtype} for subtype in SUBTYPES],
                    value=[],  # Empty by default - will analyze all subtypes
                    style={'marginBottom': '20px'}
                ),
            ], style={'maxWidth': '600px', 'margin': '0 auto', 'padding': '20px'})
        ], style={'backgroundColor': '#f8f9fa', 'borderRadius': '10px', 'padding': '20px', 'marginBottom': '30px'}),
        
        # ML Model Importance Display
        html.Div(id='gene-importance-display', 
                style={'maxWidth': '1200px', 'margin': '0 auto', 'marginBottom': '20px'}),
        
        # Results Section
        html.Div([
            html.Div(id='output-summary', style={'marginBottom': '20px'}),
            html.Div([
                dcc.Graph(id='main-plot', style={'height': '500px', 'marginBottom': '30px'}),
                dcc.Graph(id='expression-boxplot', style={'height': '400px'})
            ])
        ], style={'maxWidth': '1200px', 'margin': '0 auto', 'padding': '20px'})
    ], style={'padding': '20px', 'maxWidth': '1400px', 'margin': '0 auto'})
])


# Callback to filter gene list based on top genes checkbox
@app.callback(
    [Output('gene-selector', 'options'),
     Output('gene-selector', 'value')],
    [Input('filter-top-genes', 'value'),
     State('gene-selector', 'value')]
)
def update_gene_options(filter_enabled, current_gene):
    if filter_enabled and 'filter' in filter_enabled and TOP_GENES:
        # Filter to only top 2000 genes
        filtered_genes = [gene for gene in GENES if gene in TOP_GENES]
        options = [{'label': gene, 'value': gene} for gene in filtered_genes]
        # If current gene is not in filtered list, reset to first available
        if current_gene not in filtered_genes:
            new_value = filtered_genes[0] if filtered_genes else None
        else:
            new_value = current_gene
        return options, new_value
    else:
        # Show all genes
        options = [{'label': gene, 'value': gene} for gene in GENES]
        # If current gene is not in all genes (shouldn't happen), reset
        if current_gene not in GENES:
            new_value = GENES[0] if GENES else None
        else:
            new_value = current_gene
        return options, new_value

# Callback to display gene importance scores
@app.callback(
    Output('gene-importance-display', 'children'),
    Input('gene-selector', 'value')
)
def update_gene_importance(selected_gene):
    if not BACKEND_AVAILABLE or not selected_gene:
        return ""
    
    result = get_gene_importance(selected_gene)
    shap_score = result.get('shap_score', 'N/A')
    
    if shap_score != 'N/A':
        # Determine importance level
        if shap_score > 0.05:
            importance_level = "High"
            color = "#e74c3c"
        elif shap_score > 0.02:
            importance_level = "Medium"
            color = "#f39c12"
        else:
            importance_level = "Low"
            color = "#27ae60"
        
        return html.Div([
            html.H3("🧠 Deep Survival Neural Network Analysis", style={'textAlign': 'center', 'marginBottom': '15px', 'color': '#2c3e50'}),
            html.Div([
                html.Div([
                    html.H4(selected_gene, style={'margin': '0', 'color': '#2c3e50'}),
                    html.P(f"SHAP Importance Score: {shap_score:.5f}", 
                          style={'fontSize': '18px', 'margin': '5px 0', 'color': '#34495e'}),
                    html.P(f"Importance Level: {importance_level}", 
                          style={'fontSize': '16px', 'margin': '5px 0', 'color': color, 'fontWeight': 'bold'})
                ], style={'textAlign': 'center', 'padding': '20px', 
                         'backgroundColor': '#ffffff', 'borderRadius': '10px', 
                         'boxShadow': '0 2px 4px rgba(0,0,0,0.1)'})
            ])
        ])
    return html.Div([
        html.P(f"SHAP score not available for {selected_gene}. This gene may not be in the top 2000 genes used for model training.",
              style={'textAlign': 'center', 'color': '#7f8c8d', 'fontStyle': 'italic'})
    ])

# Callback for expression box plot
@app.callback(
    Output('expression-boxplot', 'figure'),
    Input('gene-selector', 'value')
)
def update_expression_boxplot(selected_gene):
    if not selected_gene or not BACKEND_AVAILABLE:
        return go.Figure()
    
    dataset = 'TCGA'
    result = get_expression_by_subtype(selected_gene, dataset=dataset)
    
    if 'error' in result:
        fig = go.Figure()
        fig.add_annotation(
            x=0.5, y=0.5,
            text=result['error'],
            showarrow=False,
            font=dict(size=14, color="#e67e22"),
            xref="paper", yref="paper"
        )
        fig.update_layout(
            title=f'Expression Analysis: {selected_gene}',
            xaxis=dict(showgrid=False, showticklabels=False),
            yaxis=dict(showgrid=False, showticklabels=False),
            plot_bgcolor='white',
            paper_bgcolor='white'
        )
        return fig
    
    expression_data = result['expression_data']
    
    if not expression_data:
        fig = go.Figure()
        fig.add_annotation(
            x=0.5, y=0.5,
            text="No expression data available",
            showarrow=False,
            font=dict(size=14, color="#7f8c8d"),
            xref="paper", yref="paper"
        )
        fig.update_layout(
            title=f'Expression Analysis: {selected_gene}',
            xaxis=dict(showgrid=False, showticklabels=False),
            yaxis=dict(showgrid=False, showticklabels=False),
            plot_bgcolor='white',
            paper_bgcolor='white'
        )
        return fig
    
    # Create box plot traces for each subtype
    fig = go.Figure()
    
    # Color palette for different subtypes
    colors = [
        'rgba(52, 152, 219, 0.6)',   # Blue
        'rgba(231, 76, 60, 0.6)',    # Red
        'rgba(46, 204, 113, 0.6)',   # Green
        'rgba(243, 156, 18, 0.6)',   # Orange
        'rgba(155, 89, 182, 0.6)',   # Purple
        'rgba(26, 188, 156, 0.6)'    # Teal
    ]
    line_colors = ['#3498db', '#e74c3c', '#2ecc71', '#f39c12', '#9b59b6', '#1abc9c']
    
    subtypes_sorted = sorted(expression_data.keys())
    for i, subtype in enumerate(subtypes_sorted):
        values = expression_data[subtype]
        fillcolor = colors[i % len(colors)]
        linecolor = line_colors[i % len(line_colors)]
        
        fig.add_trace(go.Box(
            y=values,
            name=subtype,
            boxmean='sd',  # Show mean and standard deviation
            marker=dict(
                color=linecolor,
                opacity=0.5,
                size=4
            ),
            line=dict(color=linecolor, width=2),
            fillcolor=fillcolor
        ))
    
    fig.update_layout(
        title=f'📊 Gene Expression Distribution by Subtype: {selected_gene}',
        xaxis_title='Subtype',
        yaxis_title='Expression Level (Z-score)',
        plot_bgcolor='white',
        paper_bgcolor='white',
        showlegend=False,
        boxmode='group',
        font=dict(size=12),
        margin=dict(l=60, r=50, t=60, b=60),
        hovermode='closest'
    )
    
    return fig

# Main analysis callback - Auto-updates when inputs change
@app.callback(
    [Output('main-plot', 'figure'),
     Output('output-summary', 'children')],
    [Input('gene-selector', 'value'),
     Input('subtype-selector', 'value')]
)
def update_analysis(selected_gene, subtypes):
    # Always use TCGA dataset
    dataset = 'TCGA'
    if not selected_gene:
        return go.Figure(), html.Div("Please select a gene to analyze.", 
                                    style={'textAlign': 'center', 'padding': '20px', 'color': '#7f8c8d'})
    
    # If no subtypes selected, analyze all available subtypes
    if not subtypes:
        subtypes = SUBTYPES if BACKEND_AVAILABLE else ['Luminal A', 'Luminal B', 'Basal-like', 'HER2-enriched', 'Normal-like']

    if not BACKEND_AVAILABLE:
        # Fallback placeholder
        fig = go.Figure()
        fig.add_trace(go.Scatter(x=[0, 12, 24], y=[1.0, 0.8, 0.6],
                                 mode='lines', name='High Expression'))
        fig.add_trace(go.Scatter(x=[0, 12, 24], y=[1.0, 0.6, 0.3],
                                 mode='lines', name='Low Expression'))
        fig.update_layout(title='Kaplan-Meier Survival Curve (Backend Not Available)',
                          xaxis_title='Months',
                          yaxis_title='Survival Probability')
        summary = html.Div([
            html.H4("Summary"),
            html.P(f"Gene: {selected_gene}"),
            html.P(f"Subtype(s): {', '.join(subtypes)}"),
            html.P("Log-rank p-value: N/A (Backend not available)")
        ])
        return fig, summary
    
    # Use real backend API - Focus on single gene, multiple subtypes
    fig = go.Figure()
    summary_items = []
    p_values = []
    
    # Plot survival curves for the selected gene across subtypes
    # Limit to 5 subtypes for clarity (can be increased if needed)
    for subtype in subtypes[:5]:
        result = get_survival_analysis(selected_gene, subtype, dataset=dataset)
        
        if 'error' in result:
            # Show user-friendly error message
            error_msg = result['error']
            if 'Dataset not loaded' in error_msg:
                error_msg = "Dataset expression data not available"
            summary_items.append(html.Div([
                html.Strong(f"{subtype}: ", style={'color': '#e74c3c'}),
                html.Span(error_msg, style={'color': '#e74c3c', 'fontSize': '14px'})
            ], style={'padding': '10px', 'marginBottom': '10px', 
                     'backgroundColor': '#fee', 'borderRadius': '5px', 'border': '1px solid #fcc'}))
            continue
        
        plot_data = result['plot_data']
        p_value = result['p_value']
        p_values.append(p_value)
        
        # Determine significance
        if p_value < 0.001:
            sig_text = "***"
            sig_color = "#e74c3c"
        elif p_value < 0.01:
            sig_text = "**"
            sig_color = "#e67e22"
        elif p_value < 0.05:
            sig_text = "*"
            sig_color = "#f39c12"
        else:
            sig_text = "ns"
            sig_color = "#95a5a6"
        
        # Add traces for high and low expression
        label_high = f"{subtype} - High Expression"
        label_low = f"{subtype} - Low Expression"
        
        fig.add_trace(go.Scatter(
            x=plot_data['high']['time'],
            y=plot_data['high']['prob'],
            mode='lines',
            name=label_high,
            line=dict(width=3, color=f"rgba({50 + len(summary_items)*30}, {100 + len(summary_items)*20}, 200, 1)" if len(summary_items) < 5 else None)
        ))
        fig.add_trace(go.Scatter(
            x=plot_data['low']['time'],
            y=plot_data['low']['prob'],
            mode='lines',
            name=label_low,
            line=dict(width=3, dash='dash', color=f"rgba({200 - len(summary_items)*30}, {50 + len(summary_items)*20}, 50, 1)" if len(summary_items) < 5 else None)
        ))
        
        summary_items.append(html.Div([
            html.Strong(f"{subtype}: ", style={'fontSize': '16px'}),
            html.Span(f"p = {p_value:.4f} ", style={'fontSize': '16px'}),
            html.Span(sig_text, style={'color': sig_color, 'fontWeight': 'bold', 'fontSize': '18px'}),
            html.Br(),
            html.Small(f"High vs Low expression comparison", style={'color': '#7f8c8d'})
        ], style={'padding': '10px', 'marginBottom': '10px', 
                 'backgroundColor': '#f8f9fa', 'borderRadius': '5px'}))
    
    fig.update_layout(
        title=f'Kaplan-Meier Survival Analysis: {selected_gene}',
        xaxis_title='Time (Months)',
        yaxis_title='Survival Probability',
        hovermode='closest',
        legend=dict(orientation="v", yanchor="top", y=1, xanchor="left", x=1.02),
        plot_bgcolor='white',
        paper_bgcolor='white',
        font=dict(size=12)
    )
    
    # Determine if analyzing all subtypes
    analyzing_all = len(subtypes) == len(SUBTYPES) if BACKEND_AVAILABLE else False
    
    summary = html.Div([
        html.H3(f"📊 Survival Analysis Results: {selected_gene}", 
                style={'textAlign': 'center', 'marginBottom': '20px', 'color': '#2c3e50'}),
        html.Div([
            html.P(f"Subtype(s): {', '.join(subtypes) if not analyzing_all else 'All Subtypes'}", 
                  style={'fontSize': '14px', 'color': '#7f8c8d', 'textAlign': 'center'}),
        ], style={'marginBottom': '20px'}),
        html.Div([
            html.P("Significance: *** p<0.001, ** p<0.01, * p<0.05, ns = not significant", 
                  style={'fontSize': '12px', 'color': '#7f8c8d', 'textAlign': 'center', 'fontStyle': 'italic'}),
        ]),
        html.Hr(),
        html.Div(summary_items)
    ])
    
    return fig, summary

# Run app
if __name__ == '__main__':
    app.run(debug=True)
