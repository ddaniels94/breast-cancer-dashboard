import dash
from dash import dcc, html, Input, Output, State
import plotly.graph_objs as go
import pandas as pd
import base64
import io

# Initialize Dash app
app = dash.Dash(__name__)
app.title = "Breast Cancer RNA Seq Explorer"

# Placeholder gene list and subtype options
GENES = ['TP53', 'BRCA1', 'BRCA2', 'ESR1', 'HER2']
SUBTYPES = ['Luminal A', 'Luminal B', 'Basal-like', 'HER2-enriched', 'Normal-like']

# Layout
app.layout = html.Div([
    html.Div([
        html.H1("🧬 Breast Cancer Survival Explorer"),
        html.Div([
            html.Label("Analysis Type", style={'marginTop': '15px'}),
            dcc.Dropdown(
                id='analysis-type',
                options=[
                    {'label': 'Survival Analysis', 'value': 'survival'},
                    {'label': 'Differential Expression', 'value': 'diffexp'},
                    {'label': 'Subtype Comparison', 'value': 'subtype'}
                ],
                value='survival',
                style={'marginBottom': '15px'}
            ),

            html.Label("Select Gene(s)", style={'marginTop': '15px'}),
            dcc.Dropdown(
                id='gene-selector',
                options=[{'label': gene, 'value': gene} for gene in GENES],
                multi=True,
                value=['TP53'],
                style={'marginBottom': '15px'}
            ),

            html.Label("Select Subtype(s)", style={'marginTop': '15px'}),
            dcc.Checklist(
                id='subtype-selector',
                options=[{'label': subtype, 'value': subtype} for subtype in SUBTYPES],
                value=['Basal-like'],
                style={'marginBottom': '15px'}
            ),

            html.Label("Expression Threshold (Percentile)"),
            dcc.Slider(
                id='threshold-slider',
                min=10, max=90, step=10, value=50,
                marks={i: f'{i}%' for i in range(10, 100, 10)}
            ),
            html.Br(),
            html.Label("Upload Custom Dataset (CSV)"),
            dcc.Upload(
                id='upload-data',
                children=html.Div(['Drag and Drop or ', html.A('Select File')]),
                style={
                    'width': '100%',
                    'height': '60px',
                    'lineHeight': '60px',
                    'borderWidth': '1px',
                    'borderStyle': 'dashed',
                    'borderRadius': '5px',
                    'textAlign': 'center',
                    'marginBottom': '20px'
                    },
        multiple=False
    ),
    dcc.Store(id='stored-data'),
            html.Button("Run Analysis", id='run-button', n_clicks=0)
        ], style={'width': '25%', 'display': 'inline-block', 'verticalAlign': 'top', 'padding': '20px'}),
    ]),

    html.Div([
        html.Div(id='output-summary'),
        dcc.Graph(id='main-plot'),
        html.Div(id='gene-table')
    ], style={'width': '70%', 'display': 'inline-block', 'padding': '20px'})
])

# Callback to parse uploaded file
@app.callback(
    Output('stored-data', 'data'),
    Input('upload-data', 'contents'),
    State('upload-data', 'filename')
)
def parse_uploaded_file(contents, filename):
    if contents is None:
        return None
    content_type, content_string = contents.split(',')
    decoded = base64.b64decode(content_string)
    try:
        if filename.endswith('.csv'):
            df = pd.read_csv(io.StringIO(decoded.decode('utf-8')))
            return df.to_dict('records')
        else:
            return None
    except Exception as e:
        print(f"Error parsing file: {e}")
        return None

# Main analysis callback
@app.callback(
    [Output('main-plot', 'figure'),
     Output('output-summary', 'children'),
     Output('gene-table', 'children')],
    [Input('run-button', 'n_clicks')],
    [State('analysis-type', 'value'),
     State('gene-selector', 'value'),
     State('subtype-selector', 'value'),
     State('threshold-slider', 'value'),
     State('stored-data', 'data')]
)
def update_analysis(n_clicks, analysis_type, genes, subtypes, threshold, stored_data):
    if n_clicks == 0:
        return go.Figure(), "", ""

    df_uploaded = pd.DataFrame(stored_data) if stored_data else None

    if analysis_type == 'survival':
        fig = go.Figure()
        fig.add_trace(go.Scatter(x=[0, 12, 24], y=[1.0, 0.8, 0.6],
                                 mode='lines', name='High Expression'))
        fig.add_trace(go.Scatter(x=[0, 12, 24], y=[1.0, 0.6, 0.3],
                                 mode='lines', name='Low Expression'))
        fig.update_layout(title='Kaplan-Meier Survival Curve',
                          xaxis_title='Months',
                          yaxis_title='Survival Probability')
        summary = html.Div([
            html.H4("Summary"),
            html.P(f"Gene(s): {', '.join(genes)}"),
            html.P(f"Subtype(s): {', '.join(subtypes)}"),
            html.P("Log-rank p-value: 0.03 (placeholder)")
        ])
        return fig, summary, ""

    elif analysis_type == 'diffexp':
        if df_uploaded is not None:
            # Example: summarize uploaded data
            df = df_uploaded.head(10)
        else:
            df = pd.DataFrame({
                'Gene': ['TP53', 'BRCA1', 'ESR1'],
                'Fold Change': [1.5, -2.1, 0.8],
                'p-value': [0.01, 0.05, 0.03]
            })
        fig = go.Figure([go.Bar(x=df['Gene'], y=df['Fold Change'])])
        fig.update_layout(title='Differential Expression (Fold Change)',
                          xaxis_title='Gene',
                          yaxis_title='Log Fold Change')
        table = html.Table([
            html.Tr([html.Th(col) for col in df.columns])] +
            [html.Tr([html.Td(df.iloc[i][col]) for col in df.columns])
             for i in range(len(df))]
        )
        return fig, "", table

    elif analysis_type == 'subtype':
        fig = go.Figure()
        fig.add_trace(go.Box(y=[2.3, 2.5, 2.7], name='Luminal A'))
        fig.add_trace(go.Box(y=[1.8, 2.0, 2.1], name='Basal-like'))
        fig.update_layout(title='Gene Expression by Subtype',
                          yaxis_title='Expression Level')
        return fig, "", ""

    return go.Figure(), "", ""

# Run app
if __name__ == '__main__':
    app.run(debug=True)
             
