{{ fullname | escape | underline }}

.. rubric:: Description

.. automodule:: {{ fullname }}

{% block classes %}
{% if classes %}
   .. rubric:: Classes

   .. autosummary::
      :toctree: .
      {% for class in classes %}
      {{ class }}
      {% endfor %}
{% endif %}
{% endblock %}

{% block functions %}
{% if functions %}
   .. rubric:: Functions

   .. autosummary::
      :toctree: .
      {% for function in functions %}
      {{ function }}
      {% endfor %}
{% endif %}
{% endblock %}