import dash
from dash import dcc, html, Input, Output, State
import plotly.graph_objs as go
import pandas as pd
 

# Initialize Dash app
app = dash.Dash(__name__)
app.title = "Breast Cancer RNA Seq Explorer"

# Placeholder gene list and subtype options
GENES = ['TP53', 'BRCA1', 'BRCA2', 'ESR1', 'HER2']
SUBTYPES = ['Luminal A', 'Luminal B', 'Basal-like', 'HER2-enriched', 'Normal-like']

# Layout
app.layout = html.Div([
    html.Div([
        html.H1("🧬 Breast Cancer RNA Seq Explorer"),
        html.Div([
            html.Label("\n","Analysis Type"),
            dcc.Dropdown(
                id='analysis-type',
                options=[
                    {'label': 'Survival Analysis', 'value': 'survival'},
                    {'label': 'Differential Expression', 'value': 'diffexp'},
                    {'label': 'Subtype Comparison', 'value': 'subtype'}
                ],
                value='survival'
            ),
            html.Label("\n","Select Gene(s)"),
            dcc.Dropdown(
                id='gene-selector',
                options=[{'label': gene, 'value': gene} for gene in GENES],
                multi=True,
                value=['TP53']
            ),
            html.Label("\n","Select Subtype(s)"),
            dcc.Checklist(
                id='subtype-selector',
                options=[{'label': subtype, 'value': subtype} for subtype in SUBTYPES],
                value=['Basal-like']
            ),
            html.Label("\n","Expression Threshold (Percentile)"),
            dcc.Slider(
                id='threshold-slider',
                min=10, max=90, step=10, value=50,
                marks={i: f'{i}%' for i in range(10, 100, 10)}
            ),
            html.Button("Run Analysis", id='run-button', n_clicks=0)
        ], style={'width': '25%', 'display': 'inline-block', 'verticalAlign': 'top', 'padding': '20px'}),
    ]),
    html.Div([
        html.Div(id='output-summary'),
        dcc.Graph(id='main-plot'),
        html.Div(id='gene-table')
    ], style={'width': '70%', 'display': 'inline-block', 'padding': '20px'})
])

# Callback
@app.callback(
    [Output('main-plot', 'figure'),
     Output('output-summary', 'children'),
     Output('gene-table', 'children')],
    [Input('run-button', 'n_clicks')],
    [State('analysis-type', 'value'),
     State('gene-selector', 'value'),
     State('subtype-selector', 'value'),
     State('threshold-slider', 'value')]
)
def update_analysis(n_clicks, analysis_type, genes, subtypes, threshold):
    if n_clicks == 0:
        return go.Figure(), "", ""

    if analysis_type == 'survival':
        # Placeholder survival analysis logic
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
        # Placeholder differential expression logic
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
        # Placeholder subtype comparison logic
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
